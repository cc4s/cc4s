/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_DEFINED
#define TCC_DEFINED

#include <tcc/Expression.hpp>
#include <tcc/Contraction.hpp>
#include <tcc/Move.hpp>
#include <tcc/Tensor.hpp>
#include <tcc/Operation.hpp>
#include <tcc/ContractionOperation.hpp>
#include <tcc/FetchOperation.hpp>
#include <tcc/IndexCounts.hpp>
#include <tcc/DryTensor.hpp>
#include <tcc/DryMachineTensor.hpp>

#include <util/StaticAssert.hpp>
#include <util/Log.hpp>

#include <vector>
#include <string>
#include <memory>

namespace tcc {
  template <typename F>
  class Tcc: public std::enable_shared_from_this<Tcc<F>> {
  protected:
    class ProtectedToken {
    };

  public:
    Tcc(
      const std::shared_ptr<MachineTensorFactory<F>> machineTensorFactory_,
      const ProtectedToken &
    ): machineTensorFactory(machineTensorFactory_) {
    }

    static std::shared_ptr<Tcc<F>> create(
      const std::shared_ptr<MachineTensorFactory<F>> machineTensorFactory_
    ) {
      return std::make_shared<Tcc<F>>(machineTensorFactory_, ProtectedToken());
    }

    /**
     * \brief Create a tcc tensor of dimensions lens_[0] x lens_[1] x ... with
     * a specified name. The underlying machine tensor will only be allocated
     * during execution of tensor operations involving this tensor.
     * Symmetryies are not supported at the moment.
     * Note that tensor objects should only be created by the Tcc object
     * which specifies the environment the tensor lives in.
     */
    std::shared_ptr<Tensor<F>> createTensor(
      const std::vector<int> &lens,
      const std::string &name
    ) {
      return std::make_shared<Tensor<F>>(
        lens, name, this->shared_from_this(),
        typename Tensor<F>::ProtectedToken()
      );
    }

    std::shared_ptr<Tensor<F>> createTensor(
      const std::shared_ptr<MachineTensor<F>> &machineTensor
    ) {
      return std::make_shared<Tensor<F>>(
        machineTensor, this->shared_from_this(),
        typename Tensor<F>::ProtectedToken()
      );
    }

    std::shared_ptr<MachineTensor<F>> createMachineTensor(
      const std::shared_ptr<Tensor<F>> &tensor
    ) {
      return machineTensorFactory->createTensor(tensor->lens, tensor->name);
    }

    std::shared_ptr<Operation<F>> compile(
      const std::shared_ptr<Expression<F>> &expression
    ) {
      auto move(std::dynamic_pointer_cast<Move<F>>(expression));
      if (move) return compile(move);
      throw new EXCEPTION("Move operation (<<=, +=, -=) expected.");
    }

    std::shared_ptr<Operation<F>> compile(
      const std::shared_ptr<Move<F>> &move
    ) {
      auto contraction(move->rhs);
      // TODO: support pure moves:
      // i.e. rhs is a contraction with exactly one element
      LOG(0, "TCC") << "compiling contraction..." << std::endl;
      LOG(2, "TCC") << "building index counts..." << std::endl;
      indexCounts = IndexCounts();
      indexCounts.add(move->lhs->indices);
      std::vector<std::shared_ptr<Operation<F>>> operations(
        contraction->factors.size()
      );
      for (unsigned int i(0); i < contraction->factors.size(); ++i) {
        operations[i] = FetchOperation<F>::create(contraction->factors[i]);
        indexCounts.add(contraction->factors[i]->indices);
      }
      triedPossibilitiesCount = 0;
      auto contractionOperation(compile(operations));

      // directly write to destination tensor instead of intermediate tensor
      contractionOperation->result = move->lhs->tensor;
      contractionOperation->resultIndices = move->lhs->indices;
      // enter the scalars alpha and beta
      contractionOperation->alpha = contraction->alpha;
      contractionOperation->beta = move->beta;

      LOG(1, "TCC") <<
        "possibilites tried=" << triedPossibilitiesCount <<
        ", FLOPS=" << contractionOperation->costs.multiplicationsCount+contractionOperation->costs.additionsCount <<
        ", maximum elements stored=" << contractionOperation->costs.maxElementsCount <<
        std::endl;
      return contractionOperation;
    }

  protected:
    /**
     * \brief Compiles the given list of Operations trying to find
     * the best order of contractions. The indexCounts are modified
     * during evaluation.
     **/
    std::shared_ptr<ContractionOperation<F>> compile(
      const std::vector<std::shared_ptr<Operation<F>>> &operations,
      const int level = 0
    ) {
      // no best contraction known at first
      std::shared_ptr<ContractionOperation<F>> bestContraction;
      for (unsigned int i(0); i < operations.size()-1; ++i) {
        std::shared_ptr<Operation<F>> a(operations[i]);
        // take out the indices of factor a
        indexCounts.add(a->getResultIndices(), -1);
        for (unsigned int j(i+1); j < operations.size(); ++j) {
          std::shared_ptr<Operation<F>> b(operations[j]);
          // take out the indices of factor b
          indexCounts.add(b->getResultIndices(), -1);

          // just compile the contraction of a&b
          std::shared_ptr<ContractionOperation<F>> abContraction(
            compile(a, b)
          );

          if (abContraction) {
            if (operations.size() == 2) {
              // we are done if there were only 2 factors to contract
              bestContraction = abContraction;
            } else {
              // otherwise, add indices of the result for further consideration
              indexCounts.add(abContraction->getResultIndices());
              // build new list of factors
              std::vector<std::shared_ptr<Operation<F>>> subOperations(
                operations.size() - 1
              );
              subOperations[0] = abContraction;
              int l(1);
              for (unsigned int k(0); k < operations.size(); ++k) {
                if (k != i && k != j) subOperations[l++] = operations[k];
              }

              // now do a recursive compilation of all the remaining factors
              std::shared_ptr<ContractionOperation<F>> fullContraction(
                compile(subOperations, level+1)
              );

              // take out indices of contraction of a&b again
              // for considering the next possibility
              indexCounts.add(abContraction->getResultIndices(), -1);

              // see if the entire contraction is currently best
              if (
                !bestContraction ||
                fullContraction->costs < bestContraction->costs
              ) {
                bestContraction = fullContraction;
                if (level == 0) { // do output only in topmost level
                  LOG(2, "TCC") <<
                    "possibilites tried=" << triedPossibilitiesCount <<
                    ", improved solution found: " <<
                    "FLOPS=" << fullContraction->costs.multiplicationsCount <<
                    ", maximum elements stored=" <<
                    fullContraction->costs.maxElementsCount << std::endl;
                }
              } else {
                if (level == 0) {
                  LOG(3, "TCC") <<
                    "possibilites tried=" << triedPossibilitiesCount <<
                    ", discarding inferior solution" << std::endl;
                }
              }
              ++triedPossibilitiesCount;
            }
          }

          // add the indices of factor b again
          indexCounts.add(b->getResultIndices());
        }
        // add the indices of factor a again
        indexCounts.add(a->getResultIndices());
      }
      return bestContraction;
    }

    std::shared_ptr<ContractionOperation<F>> compile(
      const std::shared_ptr<Operation<F>> &a,
      const std::shared_ptr<Operation<F>> &b
    ) {
      int contractedIndexDimensions[
        std::min(a->getResultIndices().length(), b->getResultIndices().length())
          + 1
      ];
      char outerIndices[
        a->getResultIndices().length() + b->getResultIndices().length() + 1
      ];
      int outerIndexDimensions[
        a->getResultIndices().length() + b->getResultIndices().length() + 1
      ];
      char uniqueIndices[
        a->getResultIndices().length() + b->getResultIndices().length() + 1
      ];
      int uniqueIndexDimensions[
        a->getResultIndices().length() + b->getResultIndices().length() + 1
      ];
      int c(0), o(0), u(0);

      // determine common unique indices
      for (unsigned int i(0); i < a->getResultIndices().length(); ++i) {
        const char index(a->getResultIndices()[i]);
        if (
          std::find(uniqueIndices, uniqueIndices+u, index) == uniqueIndices+u
        ) {
          uniqueIndices[u] = index;
          uniqueIndexDimensions[u] = a->getResult()->lens[i];
          ++u;
        }
      }
      const int uniqueAIndices(u);
      int commonIndicesCount(0);
      for (unsigned int i(0); i < b->getResultIndices().length(); ++i) {
        const char index(b->getResultIndices()[i]);
        char *previousIndex;
        if (
          (previousIndex = std::find(uniqueIndices, uniqueIndices+u, index)) ==
            uniqueIndices+u
        ) {
          uniqueIndices[u] = index;
          uniqueIndexDimensions[u] = b->getResult()->lens[i];
          ++u;
        } else if (previousIndex < uniqueIndices+uniqueAIndices) {
          ++commonIndicesCount;
        }
      }
      uniqueIndices[u] = 0;

      // skip contractions with no common indices
      if (commonIndicesCount == 0) {
        return std::shared_ptr<ContractionOperation<F>>();
      }

      int64_t outerElementsCount(1), contractedElementsCount(1);
      for (int i(0); i < u; ++i) {
        const char index(uniqueIndices[i]);
        // go through unique indices
        if (indexCounts[index] > 0) {
          // index occurs outside
          outerIndices[o] = index;
          outerElementsCount *=
            outerIndexDimensions[o] = uniqueIndexDimensions[i];
          ++o;
        } else {
          contractedElementsCount *=
            contractedIndexDimensions[c] = uniqueIndexDimensions[i];
          ++c;
        }
      }
      outerIndices[o] = 0;

      // allocate intermedate result
      std::shared_ptr<Tensor<F>> contractionResult(
        a->getResult()->getTcc()->createTensor(
          std::vector<int>(outerIndexDimensions, outerIndexDimensions+o),
          a->getResult()->getName() + b->getResult()->getName()
        )
      );
      Costs contractionCosts(
        contractionResult->getElementsCount(),
        0,
        outerElementsCount * contractedElementsCount,
        outerElementsCount * contractedElementsCount - outerElementsCount
      );
      return std::make_shared<ContractionOperation<F>>(
        a, b,
        contractionResult, static_cast<const char *>(outerIndices),
        contractionCosts,
        typename Operation<F>::ProtectedToken()
      );
    }


    std::shared_ptr<MachineTensorFactory<F>> machineTensorFactory;

    /**
     * \brief Tracks the number of occurrences of each index in the remaining
     * contraction to compile. Indices on the left hand side of the enclosing
     * assignment are also counted. This array is updated during compilation.
     **/
    IndexCounts indexCounts;

    int64_t triedPossibilitiesCount;
  };
}

#endif

