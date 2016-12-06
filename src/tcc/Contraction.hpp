/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_CONTRACTION_DEFINED
#define TCC_CONTRACTION_DEFINED

// TODO: decouple compiler from expression structure
// TODO: decouple execution from expression structure, including binding to CTF

#include <tcc/Expression.hpp>
#include <tcc/Tensor.hpp>
#include <tcc/Operation.hpp>
#include <tcc/ContractionOperation.hpp>
#include <tcc/FetchOperation.hpp>
#include <tcc/IndexCounts.hpp>
#include <util/StaticAssert.hpp>
#include <util/Log.hpp>

#include <vector>
#include <memory>

namespace tcc {
  template <typename F>
  class Contraction: public Expression<F> {
  public:
    /**
     * \brief Flattening constructor given two contractions.
     **/
    Contraction(
      const std::shared_ptr<Contraction<F>> &lhs,
      const std::shared_ptr<Contraction<F>> &rhs
    ): factors(lhs->factors) {
      factors.insert(factors.end(), rhs->factors.begin(), rhs->factors.end());
    }
    /**
     * \brief Flattening constructor given a contraction on the left hand
     * side and another expression on the right hand side.
     **/
    Contraction(
      const std::shared_ptr<Contraction<F>> &lhs,
      const std::shared_ptr<IndexedTensor<F>> &rhs
    ): factors(lhs->factors) {
      factors.push_back(rhs);
    }
    /**
     * \brief Flattening constructor given a contraction on the right hand
     * side and another expression on the left hand side.
     **/
    Contraction(
      const std::shared_ptr<IndexedTensor<F>> &lhs,
      const std::shared_ptr<Contraction<F>> &rhs
    ): factors(rhs->factors) {
      factors.push_back(lhs);
    }
    /**
     * \brief Constructor given two indexed tensors.
     **/
    Contraction(
      const std::shared_ptr<IndexedTensor<F>> &lhs,
      const std::shared_ptr<IndexedTensor<F>> &rhs
    ) {
      factors.push_back(lhs);
      factors.push_back(rhs);
    }
    /**
     * \brief Constructor given two general expressions.
     * This is currently not supported.
     **/
    Contraction(
      const std::shared_ptr<Expression<F>> &lhs,
      const std::shared_ptr<Expression<F>> &rhs
    ) {
      static_assert(
        cc4s::StaticAssert<F>::False,
        "Only contractions of contractions or tensors supported."
      );
    }
    virtual ~Contraction() {
    }

    virtual std::shared_ptr<Operation<F>> compile(
      std::string const &lhsIndices
    ) {
      LOG(0, "TCC") << "compiling contraction..." << std::endl;
      LOG(2, "TCC") << "building index counts..." << std::endl;
      indexCounts = IndexCounts();
      indexCounts.add(lhsIndices);
      std::vector<std::shared_ptr<Operation<F>>> operations(
        factors.size()
      );
      for (unsigned int i(0); i < factors.size(); ++i) {
        operations[i] = std::make_shared<FetchOperation<F>>(factors[i]);
        indexCounts.add(factors[i]->indices);
      }
      triedPossibilitiesCount = 0;
      std::shared_ptr<Operation<F>> result(compile(operations));
      LOG(1, "TCC") <<
        "possibilites tried=" << triedPossibilitiesCount <<
        ", FLOPS=" << result->costs.multiplicationsCount <<
        ", maximum elements stored=" << result->costs.maxElementsCount <<
        std::endl;
      return result;
    }

    std::vector<std::shared_ptr<IndexedTensor<F>>> factors;

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
/*
      char contractedIndices[
        std::min(a->getResultIndices().length(), b->getResultIndices().length())
          + 1
      ];
*/
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
      for (unsigned int i(0); i < b->getResultIndices().length(); ++i) {
        const char index(b->getResultIndices()[i]);
        if (
          std::find(uniqueIndices, uniqueIndices+u, index) == uniqueIndices+u
        ) {
          uniqueIndices[u] = index;
          uniqueIndexDimensions[u] = b->getResult()->lens[i];
          ++u;
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
          // index does not occur outside
//          contractedIndices[c] = index;
          contractedElementsCount *=
            contractedIndexDimensions[c] = uniqueIndexDimensions[i];
          ++c;
        }
      }
      outerIndices[o] = 0;
//      contractedIndices[c] = 0;

      // allocate intermedate result
      std::shared_ptr<Tensor<F>> contractionResult(
        std::make_shared<Tensor<F>>(
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
/*
      LOG(0, "TCC") <<
        a->getResult()->getName() << "[" << a->getResultIndices() <<
        "]*" <<
        b->getResult()->getName() << "[" << b->getResultIndices() <<
        "]: " <<
        "outerIndices=\"" << outerIndices << "\", " <<
        "contractedIndices=\"" << contractedIndices << "\"" << std::endl;
*/
      return std::make_shared<ContractionOperation<F>>(
        a, b,
        contractionResult, static_cast<const char *>(outerIndices),
        contractionCosts
      );
    }

    /**
     * \brief Tracks the number of occurrences of each index in the remaining
     * contraction to compile. Indices on the left hand side of the enclosing
     * assignment are also counted. This array is updated during compilation.
     **/
    IndexCounts indexCounts;

    int64_t triedPossibilitiesCount;
  };

  template <typename Lhs, typename Rhs>
  inline std::shared_ptr<Contraction<typename Lhs::FieldType>> operator *(
    const std::shared_ptr<Lhs> &A, const std::shared_ptr<Rhs> &B
  ) {
    static_assert(
      cc4s::TypeRelations<
        typename Lhs::FieldType, typename Rhs::FieldType
      >::Equals,
      "Only tensors of the same type can be contracted."
    );
    return std::make_shared<Contraction<typename Lhs::FieldType>>(A, B);
  }
}

#endif

