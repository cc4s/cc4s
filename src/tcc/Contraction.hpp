/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_CONTRACTION_DEFINED
#define TCC_CONTRACTION_DEFINED

#include <tcc/IndexedTensorExpression.hpp>

#include <tcc/MoveOperation.hpp>
#include <tcc/ContractionOperation.hpp>
#include <util/SharedPointer.hpp>
#include <util/StaticAssert.hpp>
#include <vector>

namespace tcc {
  template <typename F, typename TE> class Indexing;

  template <typename F, typename TE>
  class Contraction: public IndexedTensorExpression<F,TE> {
  public:
    /**
     * \brief Creates a contraction expression of the two given tensor
     * expressions A and B.
     **/
    template <typename LHS, typename RHS>
    static
    PTR(ESC(Contraction<typename RHS::FieldType, typename RHS::TensorEngine>))
    create(
      const PTR(LHS) &A, const PTR(RHS) &B
    ) {
      static_assert(
        cc4s::TypeRelations<
          typename LHS::FieldType, typename RHS::FieldType
        >::EQUALS,
        "Only tensors of the same type can be contracted."
      );
      return NEW(
        ESC(Contraction<typename RHS::FieldType, typename RHS::TensorEngine>),
        A, B,
        typename Expression<typename RHS::TensorEngine>::ProtectedToken()
      );
    }

    /**
     * \brief Creates a contraction expression of a tensor expression A
     * and a scalar alpha
     **/
    template <typename S, typename LHS>
    static
    PTR(ESC(Contraction<typename LHS::FieldType, typename LHS::TensorEngine>))
    create(
      const S alpha, const PTR(LHS) &A
    ) {
      static_assert(
        cc4s::TypeRelations<S, typename LHS::FieldType>::CASTABLE_TO,
        "The type of the scalar must be convertible to the tensor type."
      );
      return NEW(
        ESC(Contraction<typename LHS::FieldType, typename LHS::TensorEngine>),
        alpha, A,
        typename Expression<typename LHS::TensorEngine>::ProtectedToken()
      );
    }


    /**
     * \brief Flattening constructor given two contractions.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const PTR(ESC(Contraction<F,TE>)) &lhs,
      const PTR(ESC(Contraction<F,TE>)) &rhs,
      const typename Expression<TE>::ProtectedToken &
    ): alpha(lhs->alpha * rhs->alpha), factors(lhs->factors) {
      factors.insert(factors.end(), rhs->factors.begin(), rhs->factors.end());
    }
    /**
     * \brief Flattening constructor given a contraction on the left hand
     * side and another tensor result expression on the right hand side.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const PTR(ESC(Contraction<F,TE>)) &lhs,
      const PTR(ESC(IndexedTensorExpression<F,TE>)) &rhs,
      const typename Expression<TE>::ProtectedToken &
    ): alpha(lhs->alpha), factors(lhs->factors) {
      factors.push_back(rhs);
    }
    /**
     * \brief Flattening constructor given a contraction on the right hand
     * side and another tensor result expression on the left hand side.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const PTR(ESC(IndexedTensorExpression<F,TE>)) &lhs,
      const PTR(ESC(Contraction<F,TE>)) &rhs,
      const typename Expression<TE>::ProtectedToken &
    ): alpha(rhs->alpha), factors(rhs->factors) {
      factors.push_back(lhs);
    }
    /**
     * \brief Constructor given two tensor result expressions.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const PTR(ESC(IndexedTensorExpression<F,TE>)) &lhs,
      const PTR(ESC(IndexedTensorExpression<F,TE>)) &rhs,
      const typename Expression<TE>::ProtectedToken &
    ): alpha(F(1)) {
      factors.push_back(lhs);
      factors.push_back(rhs);
    }
    /**
     * \brief Constructor given a tensor result expression and a scalar.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const F alpha_,
      const PTR(ESC(IndexedTensorExpression<F,TE>)) &lhs,
      const typename Expression<TE>::ProtectedToken &
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
      const PTR(ESC(Contraction<F,TE>)) &lhs,
      const typename Expression<TE>::ProtectedToken &
    ): alpha(lhs->alpha * alpha_), factors(lhs->factors) {
    }

    virtual ~Contraction() {
    }

    virtual PTR(Operation<TE>) compile(Scope &scope) {
      std::vector<PTR(ESC(IndexedTensorOperation<F,TE>))> factorOperations(
        factors.size()
      );
      for (size_t i(0); i < factors.size(); ++i) {
        factorOperations[i] = DYNAMIC_PTR_CAST(
          ESC(IndexedTensorOperation<F,TE>), factors[i]->compile(scope)
        );
      }

      PTR(ESC(IndexedTensorOperation<F,TE>)) operation;
      if (factorOperations.size() < 2) {
        // only one operand in contraction
        operation = createSumOperation(factorOperations[0], scope);
      } else {
        // compile at least 2 contractions in best order
        operation = compileContractions(factorOperations, scope);
      }
      // enter the scaling factor alpha
      // FIXME: find proper spot for alpha
      operation->alpha = alpha;
      return operation;
    }

    // keep other overloads visible
    using Expression<TE>::compile;

    virtual void countIndices(Scope &scope) {
      for (auto &factor: factors) {
        factor->countIndices(scope);
      }
    }

    virtual operator std::string () const {
      std::stringstream stream;
      stream << "Contraction( " << alpha;
      for (auto const &factor: factors) {
        stream << ", " << std::string(*factor);
      }
      stream << " )";
      return stream.str();
    }

  protected:
    /**
     * \brief Compiles the given list of at least 2 Operations trying to
     * find the best order of contractions. The scope is modified
     * during evaluation.
     **/
    PTR(ESC(ContractionOperation<F,TE>)) compileContractions(
      const std::vector<PTR(ESC(IndexedTensorOperation<F,TE>))> &operations,
      Scope &scope,
      const unsigned int level = 0
    ) {
      // no best contraction known at first
      PTR(ESC(ContractionOperation<F,TE>)) bestContractions;
      for (unsigned int i(0); i < operations.size()-1; ++i) {
        auto a(operations[i]);
        // take out the indices of factor a
        scope.add(a->getResultIndices(), -1);
        for (unsigned int j(i+1); j < operations.size(); ++j) {
          auto b(operations[j]);
          // take out the indices of factor b
          scope.add(b->getResultIndices(), -1);

          // just compile the contraction of a&b
          auto contractionOperation(
            createContractionOperation(a, b, scope)
          );

          if (contractionOperation) {
            if (operations.size() == 2) {
              // we are done if there were only 2 factors to contract
              bestContractions = contractionOperation;
            } else {
              // otherwise, add indices of the result for further consideration
              scope.add(contractionOperation->getResultIndices());
              // build new list of factors
              std::vector<PTR(ESC(IndexedTensorOperation<F,TE>))> subOperations(
                operations.size() - 1
              );
              subOperations[0] = contractionOperation;
              int l(1);
              for (unsigned int k(0); k < operations.size(); ++k) {
                if (k != i && k != j) subOperations[l++] = operations[k];
              }

              // now do a recursive compilation of all the remaining factors
              PTR(ESC(ContractionOperation<F,TE>)) allContractions(
                compileContractions(subOperations, scope, level+1)
              );

              // take out indices of contraction of a&b again
              // for considering the next possibility
              scope.add(contractionOperation->getResultIndices(), -1);

              // see if the entire contraction is currently best
              if (
                !bestContractions ||
                allContractions->costs < bestContractions->costs
              ) {
                bestContractions = allContractions;
                if (level == 0) { // do output only in topmost level
                  LOG(2, "TCC") <<
                    "possibilites tried=" <<
                    scope.triedPossibilitiesCount <<
                    ", improved solution found: " <<
                    "FLOPS=" << allContractions->costs.multiplicationsCount <<
                    ", maximum elements stored=" <<
                    allContractions->costs.maxElementsCount << std::endl;
                }
              } else {
                if (level == 0) {
                  LOG(3, "TCC") <<
                    "possibilites tried=" <<
                    scope.triedPossibilitiesCount <<
                    ", discarding inferior solution" << std::endl;
                }
              }
              ++scope.triedPossibilitiesCount;
            }
          }

          // add the indices of factor b again
          scope.add(b->getResultIndices());
        }
        // add the indices of factor a again
        scope.add(a->getResultIndices());
      }
      return bestContractions;
    }

    /**
     * \brief Creates a ContractionOperation contracting two previously
     * compiled operations and assessing its costs.
     **/
    PTR(ESC(ContractionOperation<F,TE>)) createContractionOperation(
      const PTR(ESC(IndexedTensorOperation<F,TE>)) &a,
      const PTR(ESC(IndexedTensorOperation<F,TE>)) &b,
      Scope &scope
    ) {
      size_t contractedIndexDimensions[
        std::min(a->getResultIndices().length(), b->getResultIndices().length())
          + 1
      ];
      char outerIndices[
        a->getResultIndices().length() + b->getResultIndices().length() + 1
      ];
      size_t outerIndexDimensions[
        a->getResultIndices().length() + b->getResultIndices().length() + 1
      ];
      char uniqueIndices[
        a->getResultIndices().length() + b->getResultIndices().length() + 1
      ];
      size_t uniqueIndexDimensions[
        a->getResultIndices().length() + b->getResultIndices().length() + 1
      ];
      unsigned int c(0), o(0), u(0);

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
      const unsigned int uniqueAIndices(u);
      unsigned int commonIndicesCount(0);
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

      size_t outerElementsCount(1), contractedElementsCount(1);
      for (unsigned int i(0); i < u; ++i) {
        const char index(uniqueIndices[i]);
        // go through unique indices
        if (scope[index] > 0) {
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
      auto contractionResult(
        Tensor<F,TE>::create(
          std::vector<size_t>(outerIndexDimensions, outerIndexDimensions+o),
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
      return ContractionOperation<F,TE>::create(
        a, b,
        contractionResult, static_cast<const char *>(outerIndices),
        contractionCosts
      );
    }

    /**
     * \brief Creates a SumOperation summing
     * a compiled operation and assessing its costs.
     **/
    PTR(ESC(MoveOperation<F,TE>)) createSumOperation(
      const PTR(ESC(IndexedTensorOperation<F,TE>)) &a,
      Scope &scope
    ) {
      size_t summedIndexDimensions[a->getResultIndices().length() + 1];
      char outerIndices[a->getResultIndices().length() + 1];
      size_t outerIndexDimensions[a->getResultIndices().length() + 1];
      char uniqueIndices[a->getResultIndices().length() + 1];
      size_t uniqueIndexDimensions[a->getResultIndices().length() + 1];
      unsigned int c(0), o(0), u(0);

      // determine unique indices
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
      uniqueIndices[u] = 0;

      size_t outerElementsCount(1), summedElementsCount(1);
      for (unsigned int i(0); i < u; ++i) {
        const char index(uniqueIndices[i]);
        // go through unique indices
        if (scope[index] > 0) {
          // index occurs outside
          outerIndices[o] = index;
          outerElementsCount *=
            outerIndexDimensions[o] = uniqueIndexDimensions[i];
          ++o;
        } else {
          summedElementsCount *=
            summedIndexDimensions[c] = uniqueIndexDimensions[i];
          ++c;
        }
      }
      outerIndices[o] = 0;

      // allocate intermedate result
      auto sumResult(
        Tensor<F,TE>::create(
          std::vector<size_t>(outerIndexDimensions, outerIndexDimensions+o),
          a->getResult()->getName() + "'"
        )
      );
      // assess costs
      Costs sumCosts(
        sumResult->getElementsCount(),
        0,
        0, // no multiplications
        outerElementsCount * summedElementsCount - outerElementsCount
      );
      return MoveOperation<F,TE>::create(
        a,
        sumResult, static_cast<const char *>(outerIndices),
        sumCosts
      );
    }

    F alpha;
    std::vector<PTR(ESC(IndexedTensorExpression<F,TE>))> factors;
  };

  /**
   * \brief Creates a contraction expression of the two given tensor
   * expressions A and B using the multiplication operator *.
   **/
  template <typename LHS, typename RHS>
  inline
  PTR(ESC(Contraction<typename RHS::FieldType, typename RHS::TensorEngine>))
  operator *(
    const PTR(LHS) &A, const PTR(RHS) &B
  ) {
    return
    Contraction<typename RHS::FieldType, typename RHS::TensorEngine>::create(
      A, B
    );
  }

  /**
   * \brief Creates a contraction expression of a given tensor
   * expressions A and a scalar alpha using the multiplication operator *.
   **/
  template <typename LHS, typename S>
  inline
  PTR(ESC(Contraction<typename LHS::FieldType, typename LHS::TensorEngine>))
  operator *(
    const PTR(LHS) &A, const S alpha
  ) {
    return
    Contraction<typename LHS::FieldType, typename LHS::TensorEngine>::create(
      alpha, A
    );
  }

  /**
   * \brief Creates a contraction expression of a given tensor
   * expressions A and a scalar alpha using the multiplication operator *.
   **/
  template <typename S, typename RHS>
  inline
  PTR(ESC(Contraction<typename RHS::FieldType, typename RHS::TensorEngine>))
  operator *(
    const S alpha, const PTR(RHS) &A
  ) {
    return
    Contraction<typename RHS::FieldType, typename RHS::TensorEngine>::create(
      alpha, A
    );
  }
}

#endif

