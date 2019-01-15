/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_CONTRACTION_DEFINED
#define TCC_CONTRACTION_DEFINED

#include <tcc/Expression.hpp>

#include <tcc/ContractionOperation.hpp>
#include <tcc/MoveOperation.hpp>
#include <util/SharedPointer.hpp>
#include <util/StaticAssert.hpp>

#include <vector>

namespace tcc {
  template <typename F> class IndexedTensor;

  template <typename F>
  class Contraction: public Expression<F> {
  public:
    /**
     * \brief Creates a contraction expression of the two given tensor
     * expressions A and B.
     **/
    template <typename Lhs, typename Rhs>
    static PTR(Contraction<typename Lhs::FieldType>) create(
      const PTR(Lhs) &A, const PTR(Rhs) &B
    ) {
      static_assert(
        cc4s::TypeRelations<
          typename Lhs::FieldType, typename Rhs::FieldType
        >::EQUALS,
        "Only tensors of the same type can be contracted."
      );
      auto contraction(
        NEW(Contraction<typename Lhs::FieldType>,
          A, B,
          typename Expression<typename Lhs::FieldType>::ProtectedToken()
        )
      );
/*
      A->parent = contraction;
      B->parent = contraction;
*/
      return contraction;
    }

    /**
     * \brief Creates a contraction expression of a tensor expression A
     * and a scalar alpha
     **/
    template <typename S, typename Lhs>
    static PTR(Contraction<typename Lhs::FieldType>) create(
      const S alpha, const PTR(Lhs) &A
    ) {
      static_assert(
        cc4s::TypeRelations<S, typename Lhs::FieldType>::CASTABLE_TO,
        "The type of the scalar must be convertible to the tensor type."
      );
      auto contraction(
        NEW(Contraction<typename Lhs::FieldType>,
          alpha, A,
          typename Expression<typename Lhs::FieldType>::ProtectedToken()
        )
      );
/*
      A->parent = contraction;
*/
      return contraction;
    }


    /**
     * \brief Flattening constructor given two contractions.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const PTR(Contraction<F>) &lhs,
      const PTR(Contraction<F>) &rhs,
      const typename Expression<F>::ProtectedToken &
    ): alpha(lhs->alpha * rhs->alpha), factors(lhs->factors) {
      factors.insert(factors.end(), rhs->factors.begin(), rhs->factors.end());
    }
    /**
     * \brief Flattening constructor given a contraction on the left hand
     * side and another expression on the right hand side.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const PTR(Contraction<F>) &lhs,
      const PTR(IndexedTensor<F>) &rhs,
      const typename Expression<F>::ProtectedToken &
    ): alpha(lhs->alpha), factors(lhs->factors) {
      factors.push_back(rhs);
    }
    /**
     * \brief Flattening constructor given a contraction on the right hand
     * side and another expression on the left hand side.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const PTR(IndexedTensor<F>) &lhs,
      const PTR(Contraction<F>) &rhs,
      const typename Expression<F>::ProtectedToken &
    ): alpha(rhs->alpha), factors(rhs->factors) {
      factors.push_back(lhs);
    }
    /**
     * \brief Constructor given two indexed tensors.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const PTR(IndexedTensor<F>) &lhs,
      const PTR(IndexedTensor<F>) &rhs,
      const typename Expression<F>::ProtectedToken &
    ): alpha(F(1)) {
      factors.push_back(lhs);
      factors.push_back(rhs);
    }
    /**
     * \brief Constructor given two general expressions.
     * This is currently not supported.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const PTR(Expression<F>) &lhs,
      const PTR(Expression<F>) &rhs,
      const typename Expression<F>::ProtectedToken &
    ) {
      static_assert(
        cc4s::StaticAssert<F>::FALSE,
        "Only contractions of contractions or tensors supported."
      );
    }
    /**
     * \brief Constructor given an indexed tensor and a scalar.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const F alpha_,
      const PTR(IndexedTensor<F>) &lhs,
      const typename Expression<F>::ProtectedToken &
    ): alpha(alpha_) {
      factors.push_back(lhs);
    }
    /**
     * \brief Falttening constructor given a contraction and a scalar.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const F alpha_,
      const PTR(Contraction<F>) &lhs,
      const typename Expression<F>::ProtectedToken &
    ): alpha(lhs->alpha * alpha_), factors(lhs->factors) {
    }

    virtual ~Contraction() {
    }

    virtual PTR(Operation<F>) compile(IndexCounts &indexCounts) {
      std::vector<PTR(Operation<F>)> factorOperations(factors.size());
      for (unsigned int i(0); i < factors.size(); ++i) {
        factorOperations[i] = factors[i]->compile(indexCounts);
      }

      PTR(TensorResultOperation<F>) operation;
      if (factorOperations.size() < 2) {
        // only one operand in contraction: do move directly
        operation = createMoveOperation(factorOperations[0]);
      } else {
        // compile at least 2 contractions in correct order
        operation = compileContractions(factorOperations, indexCounts);
      }
      // enter the scaling factor alpha
      operation->alpha = alpha;
      return operation;
    }

    virtual void countIndices(IndexCounts &indexCounts) {
      for (auto &factor: factors) {
        factor->countIndices(indexCounts);
      }
    }

  protected:
    /**
     * \brief Compiles the given list of at least 2 Operations trying to
     * find the best order of contractions. The indexCounts are modified
     * during evaluation.
     **/
    PTR(ContractionOperation<F>) compileContractions(
      const std::vector<PTR(Operation<F>)> &operations,
      IndexCounts &indexCounts,
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
            createContractionOperation(a, b, indexCounts)
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
                compileContractions(subOperations, indexCounts, level+1)
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
                    "possibilites tried=" <<
                    indexCounts.triedPossibilitiesCount <<
                    ", improved solution found: " <<
                    "FLOPS=" << allContractions->costs.multiplicationsCount <<
                    ", maximum elements stored=" <<
                    allContractions->costs.maxElementsCount << std::endl;
                }
              } else {
                if (level == 0) {
                  LOG(3, "TCC") <<
                    "possibilites tried=" <<
                    indexCounts.triedPossibilitiesCount <<
                    ", discarding inferior solution" << std::endl;
                }
              }
              ++indexCounts.triedPossibilitiesCount;
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
    PTR(MoveOperation<F>) createMoveOperation(
      const PTR(Operation<F>) &a
    ) {
      // allocate intermedate result assuming identical indices as argument
      PTR(Tensor<F>) moveResult(
        a->getResult()->getTcc()->createTensor(
          a->getResult(), a->getResult()->getName() + "'"
        )
      );
      int64_t elementsCount(1);
      for (int i: a->getResult()->lens) {
        elementsCount *= i;
      }
      // assess costs
      Costs moveCosts(
        moveResult->getElementsCount(),
        moveResult->getElementsCount(),
        0,  // no multiplications
        moveResult->getElementsCount() // number of additions
      );
      return MoveOperation<F>::create(
        a,
        moveResult, a->getResultIndices().c_str(),
        moveCosts
      );
    }

    /**
     * \brief Creates a ContractionOperation contracting two previously
     * compiled operations and assessing its costs.
     **/
    PTR(ContractionOperation<F>) createContractionOperation(
      const PTR(Operation<F>) &a,
      const PTR(Operation<F>) &b,
      IndexCounts &indexCounts
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

    F alpha;
    std::vector<PTR(IndexedTensor<F>)> factors;

    friend class Tcc<F>;
  };

  /**
   * \brief Creates a contraction expression of the two given tensor
   * expressions A and B using the multiplication operator *.
   **/
  template <typename Lhs, typename Rhs>
  inline PTR(Contraction<typename Lhs::FieldType>) operator *(
    const PTR(Lhs) &A, const PTR(Rhs) &B
  ) {
    return Contraction<typename Lhs::FieldType>::create(A, B);
  }

  /**
   * \brief Creates a contraction expression of a given tensor
   * expressions A and a scalar alpha using the multiplication operator *.
   **/
  template <typename Lhs, typename S>
  inline PTR(Contraction<typename Lhs::FieldType>) operator *(
    const PTR(Lhs) &A, const S alpha
  ) {
    return Contraction<typename Lhs::FieldType>::create(alpha, A);
  }

  /**
   * \brief Creates a contraction expression of a given tensor
   * expressions A and a scalar alpha using the multiplication operator *.
   **/
  template <typename S, typename Rhs>
  inline PTR(Contraction<typename Rhs::FieldType>) operator *(
    const S alpha, const PTR(Rhs) &A
  ) {
    return Contraction<typename Rhs::FieldType>::create(alpha, A);
  }
}

#endif

