/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_DEFINED
#define TCC_DEFINED

#include <tcc/IndexedTensor.hpp>
#include <tcc/Move.hpp>
#include <tcc/Contraction.hpp>
#include <tcc/Sequence.hpp>
#include <tcc/Tensor.hpp>
#include <tcc/FetchOperation.hpp>
#include <tcc/MoveOperation.hpp>
#include <tcc/ContractionOperation.hpp>
#include <tcc/OperationSequence.hpp>
#include <tcc/IndexCounts.hpp>

#include <util/SharedPointer.hpp>
#include <util/Log.hpp>

#include <vector>
#include <string>
#include <memory>

// FIXME: tensor <<= sclar * tensor throws static exception
// TODO: *= operator with scalars
// TODO: heuristics: limit number of simultaneously considered intermediates
// TODO: fix max memory assessment
// TODO: mixed type tensor operations
// TODO: unary and binary function application
// TODO: expression definitions with local index renaming
// TODO: permutation and anti-permutation operator
// TODO: common subexpression optimization
// TODO: support hard memory limit for costs
// TODO: support slicing and looping over indices for memory reduction

namespace tcc {
  template <typename F>
  class Tcc: public std::enable_shared_from_this<Tcc<F>> {
  protected:
    class ProtectedToken {
    };

  public:
    Tcc(
      const PTR(MachineTensorFactory<F>) machineTensorFactory_,
      const ProtectedToken &
    ): machineTensorFactory(machineTensorFactory_) {
    }

    static PTR(Tcc<F>) create(
      const PTR(MachineTensorFactory<F>) machineTensorFactory_
    ) {
      return NEW(Tcc<F>, machineTensorFactory_, ProtectedToken());
    }

    /**
     * \brief Create a tcc tensor of dimensions lens_[0] x lens_[1] x ... with
     * a specified name. The underlying machine tensor will only be allocated
     * during execution of tensor operations involving this tensor.
     * Symmetryies are not supported at the moment.
     * Note that tensor objects should only be created by the Tcc object
     * which specifies the environment the tensor lives in.
     **/
    PTR(Tensor<F>) createTensor(
      const std::vector<int> &lens,
      const std::string &name
    ) {
      return NEW(Tensor<F>,
        lens, name, this->shared_from_this(),
        typename Tensor<F>::ProtectedToken()
      );
    }

    /**
     * \brief Create an empty tensor of identical shape as the given tensor.
     * The name, however, should differ.
     **/
    PTR(Tensor<F>) createTensor(
      const PTR(Tensor<F>) &tensor,
      const std::string &name
    ) {
      return createTensor(tensor->lens, name);
    }

    /**
     * \brief Create a tcc tensor that uses a given machine tensor as its
     * machine tensor for representing the tensor data.
     **/
    PTR(Tensor<F>) createTensor(
      const PTR(MachineTensor<F>) &machineTensor
    ) {
      return NEW(Tensor<F>,
        machineTensor, this->shared_from_this(),
        typename Tensor<F>::ProtectedToken()
      );
    }

    PTR(MachineTensor<F>) createMachineTensor(
      const PTR(Tensor<F>) &tensor
    ) {
      return machineTensorFactory->createTensor(tensor->lens, tensor->name);
    }

    PTR(Operation<F>) compile(
      const PTR(Expression<F>) &expression
    ) {
      auto sequence(dynamic_cast<PTR(Sequence<F>)>(expression));
      if (sequence) return compile(sequence);
      auto move(dynamic_cast<PTR(Move<F>)>(expression));
      if (move) return compile(move);
      throw new EXCEPTION("Sequence (,) of move operation (<<=, +=, -=) expected.");
    }

    PTR(Sequence<F>) emptySequence() {
      return NEW(Sequence<F>);
    }

    PTR(Operation<F>) compile(
      const PTR(Sequence<F>) &sequence
    ) {
      std::vector<PTR(Operation<F>)> operations(
        sequence->moves.size()
      );
      for (unsigned int i(0); i < sequence->moves.size(); ++i) {
        operations[i] = compile(sequence->moves[i]);
      }
      return OperationSequence<F>::create(operations);
    }

    PTR(Operation<F>) compile(
      const PTR(Move<F>) &move
    ) {
      LOG(2, "TCC") << "compiling contraction..." << std::endl;

      indexCounts = IndexCounts();
      indexCounts.add(move->lhs->indices);
      auto contraction(move->rhs);
      std::vector<PTR(Operation<F>)> operations(
        contraction->factors.size()
      );
      for (unsigned int i(0); i < contraction->factors.size(); ++i) {
        operations[i] = FetchOperation<F>::create(contraction->factors[i]);
        indexCounts.add(contraction->factors[i]->indices);
      }
      triedPossibilitiesCount = 0;

      PTR(TensorResultOperation<F>) operation;
      if (operations.size() < 2) {
        operation = MoveOperation<F>::create(
          operations[0],
          move->lhs->tensor, move->lhs->indices.c_str(),
          Costs(move->lhs->tensor->getElementsCount())
        );
      } else {
        operation = compileContractions(operations);
        // directly write to destination tensor instead of intermediate tensor
        operation->result = move->lhs->tensor;
        operation->resultIndices = move->lhs->indices;
      }
      // enter the scalars alpha and beta
      operation->alpha = contraction->alpha;
      operation->beta = move->beta;

      LOG(2, "TCC") <<
        "possibilites tried=" << triedPossibilitiesCount <<
        ", FLOPS=" <<
          operation->costs.multiplicationsCount +
          operation->costs.additionsCount <<
        ", maximum elements stored=" <<
          operation->costs.maxElementsCount <<
        std::endl;

      return operation;
    }

  protected:
    /**
     * \brief Compiles the given list of at least 2 Operations trying to
     * find the best order of contractions. The indexCounts are modified
     * during evaluation.
     **/
    PTR(ContractionOperation<F>) compileContractions(
      const std::vector<PTR(Operation<F>)> &operations,
      const int level = 0
    ) {
      // no best contraction known at first
      PTR(ContractionOperation<F>) bestContractions;
      for (unsigned int i(0); i < operations.size()-1; ++i) {
        PTR(Operation<F>) a(operations[i]);
        // take out the indices of factor a
        indexCounts.add(a->getResultIndices(), -1);
        for (unsigned int j(i+1); j < operations.size(); ++j) {
          PTR(Operation<F>) b(operations[j]);
          // take out the indices of factor b
          indexCounts.add(b->getResultIndices(), -1);

          // just compile the contraction of a&b
          PTR(ContractionOperation<F>) contractionOperation(
            createContractionOperation(a, b)
          );

          if (contractionOperation) {
            if (operations.size() == 2) {
              // we are done if there were only 2 factors to contract
              bestContractions = contractionOperation;
            } else {
              // otherwise, add indices of the result for further consideration
              indexCounts.add(contractionOperation->getResultIndices());
              // build new list of factors
              std::vector<PTR(Operation<F>)> subOperations(
                operations.size() - 1
              );
              subOperations[0] = contractionOperation;
              int l(1);
              for (unsigned int k(0); k < operations.size(); ++k) {
                if (k != i && k != j) subOperations[l++] = operations[k];
              }

              // now do a recursive compilation of all the remaining factors
              PTR(ContractionOperation<F>) allContractions(
                compileContractions(subOperations, level+1)
              );

              // take out indices of contraction of a&b again
              // for considering the next possibility
              indexCounts.add(contractionOperation->getResultIndices(), -1);

              // see if the entire contraction is currently best
              if (
                !bestContractions ||
                allContractions->costs < bestContractions->costs
              ) {
                bestContractions = allContractions;
                if (level == 0) { // do output only in topmost level
                  LOG(2, "TCC") <<
                    "possibilites tried=" << triedPossibilitiesCount <<
                    ", improved solution found: " <<
                    "FLOPS=" << allContractions->costs.multiplicationsCount <<
                    ", maximum elements stored=" <<
                    allContractions->costs.maxElementsCount << std::endl;
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
      return bestContractions;
    }

    /**
     * \brief Creates a ContractionOperation contracting two previously
     * compiled operations and assessing its costs.
     **/
    PTR(ContractionOperation<F>) createContractionOperation(
      const PTR(Operation<F>) &a,
      const PTR(Operation<F>) &b
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
      PTR(Tensor<F>) contractionResult(
        a->getResult()->getTcc()->createTensor(
          std::vector<int>(outerIndexDimensions, outerIndexDimensions+o),
          a->getResult()->getName() + b->getResult()->getName()
        )
      );
      // assess costs
      Costs contractionCosts(
        contractionResult->getElementsCount(),
        0,
        outerElementsCount * contractedElementsCount,
        outerElementsCount * contractedElementsCount - outerElementsCount
      );
      return ContractionOperation<F>::create(
        a, b,
        contractionResult, static_cast<const char *>(outerIndices),
        contractionCosts
      );
    }


    PTR(MachineTensorFactory<F>) machineTensorFactory;

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

