/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_CONTRACTION_DEFINED
#define TENSOR_CONTRACTION_DEFINED

#include <tcc/TensorExpression.hpp>
#include <tcc/DryTensor.hpp>
#include <tcc/TensorOperation.hpp>
#include <tcc/TensorContractionOperation.hpp>
#include <tcc/TensorFetchOperation.hpp>
#include <tcc/IndexCounts.hpp>
#include <util/StaticAssert.hpp>
#include <util/Log.hpp>
#include <vector>

namespace cc4s {
  template <typename F>
  class TensorContraction: public TensorExpression<F> {
  public:
    /**
     * \brief Flattening constructor given two contractions.
     **/
    TensorContraction(
      TensorContraction<F> *lhs, TensorContraction<F> *rhs
    ): factors(lhs.factors) {
      factors.insert(factors.end(), rhs->factors.begin(), rhs->factors.end());
      // the factors from both lhs and rhs contractions are now contained here
      lhs->factors.clear(); rhs->factors.clear();
      delete lhs, rhs;
    }
    /**
     * \brief Flattening constructor given a contraction on the left hand
     * side and another expression on the right hand side.
     **/
    TensorContraction(
      TensorContraction<F> *lhs, IndexedTensor<F> *rhs
    ): factors(lhs->factors) {
      factors.push_back(rhs);
      // the factors from the lhs expression are now contained here
      lhs->factors.clear();
      delete lhs;
    }
    /**
     * \brief Flattening constructor given a contraction on the right hand
     * side and another expression on the left hand side.
     **/
    TensorContraction(
      IndexedTensor<F> *lhs, TensorContraction<F> *rhs
    ): factors(rhs.factors) {
      factors.push_back(lhs);
      // the factors from the rhs expression are now contained here
      rhs->factors.clear();
      delete rhs;
    }
    /**
     * \brief Constructor given two indexed tensors.
     **/
    TensorContraction(
      IndexedTensor<F> *lhs, IndexedTensor<F> *rhs
    ) {
      factors.push_back(lhs);
      factors.push_back(rhs);
    }
    /**
     * \brief Constructor given two general expressions.
     * This is currently not supported.
     **/
    TensorContraction(
      TensorExpression<F> *lhs, TensorExpression<F> *rhs
    ) {
      static_assert(
        StaticAssert<F>::False,
        "Only contractions of contractions or tensors supported."
      );
    }
    virtual ~TensorContraction() {
      // subexpressions are dependent entities: delete each factor
      for (auto factor(factors.begin()); factor != factors.end(); ++factor) {
        delete *factor;
      }
    }

    virtual void log() const {
      for (auto factor(factors.begin()); factor != factors.end(); ++factor) {
        (*factor)->log();
      }
      LOG(0, "TCC") << factors.size() << " tensors contracted" << std::endl;
    }

    virtual TensorOperation<F> *compile(std::string const &lhsIndices) {
      LOG(0, "TCC") << "compiling contraction..." << std::endl;
      LOG(0, "TCC") << "building index counts..." << std::endl;
      indexCounts = IndexCounts();
      indexCounts.add(lhsIndices);
      std::vector<TensorOperation<F> *> operations(factors.size());
      for (int i(0); i < factors.size(); ++i) {
        operations[i] = new TensorFetchOperation<F>(factors[i]);
        indexCounts.add(factors[i]->indices);
      }
      return compile(operations);
    }

    std::vector<IndexedTensor<F> *> factors;

  protected:
    /**
     * \brief Compiles the given list of TensorOperations trying to find
     * the best order of contractions. The indexCounts are modified
     * during evaluation.
     **/
    TensorContractionOperation<F> *compile(
      std::vector<TensorOperation<F> *> operations
    ) {
      // no best contraction known at first
      TensorContractionOperation<F> *bestContraction(nullptr);
      for (int i(0); i < operations.size()-1; ++i) {
        TensorOperation<F> *a(operations[i]);
        // take out the indices of factor a
        indexCounts.add(a->getResultIndices(), -1);
        for (int j(i+1); j < operations.size(); ++j) {
          TensorOperation<F> *b(operations[j]);
          // take out the indices of factor b
          indexCounts.add(b->getResultIndices(), -1);

          // just compile the contraction of a&b
          TensorContractionOperation<F> *abContraction(compile(a, b));

          if (operations.size() == 2) {
            // we are done if there were only 2 factors to contract
            bestContraction = abContraction;
          } else {
            // otherwise, add indices of the result for further consideration
            indexCounts.add(abContraction->getResultIndices());
            // build new list of factors
            std::vector<TensorOperation<F> *> subOperations(
              operations.size() - 1
            );
            subOperations[0] = abContraction;
            int l(1);
            for (int k(0); k < operations.size(); ++k) {
              if (k != i && k != j) subOperations[l++] = operations[k];
            }

            // now do a recursive compilation of all the remaining factors
            TensorContractionOperation<F> *fullContraction(
              compile(subOperations)
            );

            // take out indices of contraction of a&b again
            // for considering the next possibility
            indexCounts.add(abContraction->getResultIndices(), -1);

            // see if the entire contraction is currently best
            if (
              !bestContraction ||
              fullContraction->costs < bestContraction->costs
            ) {
              LOG(0, "TCC") << "improved solution found: " <<
                fullContraction->costs.multiplicationsCount << " FLOPS, " <<
                fullContraction->costs.maxElementsCount << " elements" <<
                std::endl;
              if (bestContraction) {
                // yes: throw away the previous best but keep the fetch ops
                if (bestContraction->clearLeavingFetches()) {
                  delete bestContraction;
                }
              }
              bestContraction = fullContraction;
            } else {
              LOG(0, "TCC") << "discarding inferior solution" << std::endl;
              // no: throw away this one but keep the fetch operations
              if (fullContraction->clearLeavingFetches()) {
                delete fullContraction;
              }
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

    TensorContractionOperation<F> *compile(
      TensorOperation<F> *a, TensorOperation<F> *b
    ) {
      char contractedIndices[
        std::min(a->getResultIndices().length(), b->getResultIndices().length())
          + 1
      ];
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
      int outerIndexSymmetries[
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
      for (int i(0); i < a->getResultIndices().length(); ++i) {
        const char index(a->getResultIndices()[i]);
        if (
          std::find(uniqueIndices, uniqueIndices+u, index) == uniqueIndices+u
        ) {
          uniqueIndices[u] = index;
          uniqueIndexDimensions[u] = a->getResult()->lens[i];
          ++u;
        }
      }
      for (int i(0); i < b->getResultIndices().length(); ++i) {
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
          outerIndexSymmetries[o] = 0;
          ++o;
        } else {
          // index does not occur outside
          contractedIndices[c] = index;
          contractedElementsCount *=
            contractedIndexDimensions[c] = uniqueIndexDimensions[i];
          ++c;
        }
      }
      outerIndices[o] = 0;
      contractedIndices[c] = 0;

      // allocate intermedate result
      DryTensor<F> *contractionResult(
        new DryTensor<F>(o, outerIndexDimensions, outerIndexSymmetries)
      );
      contractionResult->set_name(
        a->getResult()->get_name() + b->getResult()->get_name()
      );
      // TODO: name intermediate result tensor
      Costs contractionCosts(
        contractionResult->getElementsCount(),
        0,
        outerElementsCount * contractedElementsCount,
        outerElementsCount * contractedElementsCount - outerElementsCount
      );

      LOG(0, "TCC") <<
        a->getResult()->get_name() << "[" << a->getResultIndices() <<
        "]*" <<
        b->getResult()->get_name() << "[" << b->getResultIndices() <<
        "]: " <<
        "outerIndices=\"" << outerIndices << "\", " <<
        "contractedIndices=\"" << contractedIndices << "\"" << std::endl;

      return new TensorContractionOperation<F>(
        a, b,
        contractionResult, outerIndices,
        contractionCosts
      );
    }

    /**
     * \brief Tracks the number of occurrences of each index in the remaining
     * contraction to compile. Indices on the left hand side of the enclosing
     * assignment are also counted. This array is updated during compilation.
     **/
    IndexCounts indexCounts;
  };

  template <typename Lhs, typename Rhs>
  inline TensorContraction<typename Lhs::FieldType> &operator *(
    Lhs &A, Rhs &B
  ) {
    static_assert(
      TypeRelations<typename Lhs::FieldType, typename Rhs::FieldType>::Equals,
      "Only tensors of the same type can be contracted."
    );
    return *new TensorContraction<typename Lhs::FieldType>(&A, &B);
  }
}

#endif

