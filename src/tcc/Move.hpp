/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_MOVE_DEFINED
#define TCC_MOVE_DEFINED

#include <tcc/IndexedTensorExpression.hpp>

#include <tcc/Contraction.hpp>
#include <tcc/MoveOperation.hpp>

#include <util/SharedPointer.hpp>
#include <util/StaticAssert.hpp>

namespace tcc {
  template <typename F, typename TE>
  class Indexing;

  template <typename F, typename TE>
  class Move: public IndexedTensorExpression<F,TE> {
  public:
    /**
     * \brief Moves the given right hand side expression into the given
     * indexed tensor, returning the indexed tensor as result expression
     * for possible further operations.
     **/
    template <typename LHS, typename RHS, typename S>
    static inline
    PTR(ESC(Move<typename LHS::FieldType, typename LHS::TensorEngine>)) create(
      const PTR(LHS) &lhs,
      const PTR(RHS) &rhs,
      const S beta
    ) {
      static_assert(
        cc4s::TypeRelations<S, typename LHS::FieldType>::CASTABLE_TO,
        "The type of the scalar must be convertible to the tensor type."
      );
      return NEW(ESC(Move<typename LHS::FieldType,TE>),
        lhs, rhs, beta,
        typename Expression<TE>::ProtectedToken()
      );
    }

    /**
     * \brief Creates a move expression of the right hand side tensor
     * expression rhs in the left hand tensor expression lhs.
     * Not indended for direct invocation. Use Move::create instead
     **/
/*
    FIXME: find exclusion
    Move(
      const PTR(ESC(IndexedTensor<F,TE>)) &lhs_,
      const PTR(Expression<TE>) &rhs_,
      const F beta_,
      const typename Expression<TE>::ProtectedToken &
    ) {
      static_assert(
        cc4s::StaticAssert<F>::FALSE,
        "The result type of the right-hand-side must match that of the left-hand-side."
      );
    }
*/

    /**
     * \brief Creates a move expression where the right hand side is
     * a indexed tensor expression but not a contraction.
     * It will be wrapped in a Contraction
     * Expression so that all moves have contractions as their right hand side.
     * Not indended for direct invocation. Use Move::create instead
     **/
    Move(
      const PTR(ESC(Indexing<F,TE>)) &lhs_,
      const PTR(ESC(IndexedTensorExpression<F,TE>)) &rhs_,
      const F beta_,
      const typename Expression<TE>::ProtectedToken &
    ): lhs(lhs_), rhs(Contraction<F,TE>::create(1, rhs_)), beta(beta_) {
    }

    /**
     * \brief Creates a move expression where the right hand side is
     * a Contraction expression.
     * Not indended for direct invocation. Use Move::create instead
     **/
    Move(
      const PTR(ESC(Indexing<F,TE>)) &lhs_,
      const PTR(ESC(Contraction<F,TE>)) &rhs_,
      const F beta_,
      const typename Expression<TE>::ProtectedToken &
    ): lhs(lhs_), rhs(rhs_), beta(beta_) {
    }
    virtual ~Move() {
    }

    // each move has its private index namespace so disregard the outer
    // indexCounts
    virtual PTR(Operation<TE>) compile(IndexCounts &) {
      LOG(2, "TCC") << "compiling move..." << std::endl;
      // create a new namespace of indices
      IndexCounts indexCounts;
      // and determine how often each index is used
      countIndices(indexCounts);

      // TODO: currently only move(lhs,contraction(factors...)) are done
      indexCounts.triedPossibilitiesCount = 0;
      auto operation(
        // compile right-hand-side in the namespace of this move
        DYNAMIC_PTR_CAST(
          ESC(IndexedTensorOperation<F,TE>), rhs->compile(indexCounts)
        )
      );

      assumeOrCheckShape(operation);

      // write operation results directly to lhs tensor instead of intermediate
      operation->result = lhs->source;
      operation->resultIndices = lhs->indices;
      // enter the beta factor of this move
      // TODO: enter in contraction or move
      operation->beta = beta;

      LOG(2, "TCC") <<
        "possibilites tried=" << indexCounts.triedPossibilitiesCount <<
        ", FLOPS=" <<
          operation->costs.multiplicationsCount +
          operation->costs.additionsCount <<
        ", maximum elements stored=" <<
          operation->costs.maxElementsCount <<
        std::endl;

      return operation;
    }

    virtual void countIndices(IndexCounts &indexCounts) {
      lhs->countIndices(indexCounts);
      rhs->countIndices(indexCounts);
    }

  protected:
    void assumeOrCheckShape(
      const PTR(ESC(IndexedTensorOperation<F,TE>)) &operation
    ) {
      if (lhs->indices.length() != operation->getResultIndices().length()) {
        throw new EXCEPTION(
          "Number of indices of left-hand-side tensor " + lhs->source->getName() +
          " must match the number of indices of the result tensor " +
          operation->getResult()->getName()
        );
      }
      // determine shape of lhs
      std::vector<size_t> lens(operation->getResult()->getLens().size());
      size_t lenOfIndex[std::numeric_limits<char>::max()+1];
      // enter which indices correspond to which length on the rhs
      for (unsigned int i(0); i < lhs->indices.length(); ++i) {
        lenOfIndex[static_cast<unsigned int>(operation->getResultIndices()[i])]=
          operation->result->getLens()[i];
      }
      // use for lhs
      for (unsigned int i(0); i < lhs->indices.length(); ++i) {
        lens[i] = lenOfIndex[static_cast<unsigned int>(lhs->indices[i])];
      }
      if (!lhs->source->assumedShape) {
        // assume
        lhs->source->lens = lens;
        lhs->source->assumedShape = true;
        // or check shape
      } else if (lhs->source->getLens() != lens) {
        if (lhs->source->lens != operation->getResult()->getLens()) {
          throw new EXCEPTION(
            "Shape of left-hand-side tensor " + lhs->source->getName() +
            " must match the shape of the result tensor " +
            operation->getResult()->getName()
          );
        }
      }
    }

    // TODO: do lvalueCompile for general left-hand-side compilation
    PTR(ESC(Indexing<F,TE>)) lhs;
    PTR(ESC(Contraction<F,TE>)) rhs;
    F beta;
  };

  /**
   * \brief Creates an updating move expression.
   **/
  template <typename RHS>
  inline
  PTR(ESC(Move<typename RHS::FieldType,typename RHS::TensorEngine>))
  operator +=(
    const
    PTR(
      ESC(Indexing<typename RHS::FieldType,typename RHS::TensorEngine>)
    ) &lhs,
    const PTR(RHS) &rhs
  ) {
    return Move<typename RHS::FieldType,typename RHS::TensorEngine>::create(
      lhs, rhs, typename RHS::FieldType(1)
    );
  }
  /**
   * \brief Creates an updating move expression.
   **/
  template <typename RHS>
  inline
  PTR(ESC(Move<typename RHS::FieldType,typename RHS::TensorEngine>))
  operator -=(
    const
    PTR(
      ESC(Indexing<typename RHS::FieldType,typename RHS::TensorEngine>)
    ) &lhs,
    const PTR(RHS) &rhs
  ) {
    return Move<typename RHS::FieldType,typename RHS::TensorEngine>::create(
      lhs,
      Contraction<typename RHS::FieldType,typename RHS::TensorEngine>::create(
        typename RHS::FieldType(-1), rhs
      ),
      typename RHS::FieldType(1)
    );
  }

  /**
   * \brief Moves the given right hand side expression into the given
   * indexed tensor, returning the indexed tensor as result expression
   * for possible further operations.
   * Note that the operator = cannot be used since all expressions are
   * represented by shared pointers.
   **/
  template <typename RHS>
  inline
  PTR(ESC(Move<typename RHS::FieldType, typename RHS::TensorEngine>))
  operator <<=(
    const
    PTR(
      ESC(Indexing<typename RHS::FieldType,typename RHS::TensorEngine>)
    ) &lhs,
    const PTR(RHS) &rhs
  ) {
    return Move<typename RHS::FieldType, typename RHS::TensorEngine>::create(
      lhs, rhs, typename RHS::FieldType(0)
    );
  }

  template <typename LHS, typename RHS>
  inline
  PTR(ESC(Move<typename RHS::FieldType, typename RHS::TensorEngine>))
  operator <<=(
    const PTR(LHS) &, const PTR(RHS) &rhs
  ) {
    static_assert(
      cc4s::TypeRelations<
        typename LHS::FieldType, typename RHS::FieldType
      >::EQUALS,
      "Move operations requires tensors of same type."
    );
    static_assert(
      cc4s::StaticAssert<RHS>::FALSE,
      "Only indexed tensors may be used as the left hand side of a move operation."
    );
    return PTR(ESC(Move<typename RHS::FieldType,typename RHS::TensorEngine>))();
  }

  template <typename F, typename TE>
  inline std::ostream &operator <<(
    std::ostream &stream, const Indexing<F,TE> &t
  ) {
    return stream << t.getTensor()->getName() <<
      "[" << t->getResultIndices() << "]";
  }
}

#endif

