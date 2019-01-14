/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_MOVE_DEFINED
#define TCC_MOVE_DEFINED

#include <tcc/Expression.hpp>
#include <tcc/Contraction.hpp>

#include <util/SharedPointer.hpp>
#include <util/StaticAssert.hpp>

namespace tcc {
  template <typename F>
  class IndexedTensor;

  template <typename F>
  class Move: public Expression<F> {
  public:
    /**
     * \brief Moves the given right hand side expression into the given
     * indexed tensor, returning the indexed tensor as result expression
     * for possible further operations.
     **/
    template <typename Lhs, typename Rhs, typename S>
    static inline PTR(Move<typename Lhs::FieldType>) create(
      const PTR(Lhs) &lhs,
      const PTR(Rhs) &rhs,
      const S beta
    ) {
      static_assert(
        cc4s::TypeRelations<S, typename Lhs::FieldType>::CASTABLE_TO,
        "The type of the scalar must be convertible to the tensor type."
      );
      auto move(
        NEW(Move<typename Lhs::FieldType>,
          lhs, rhs, beta,
          typename Expression<typename Lhs::FieldType>::ProtectedToken()
        )
      );
      lhs->parent = move;
      rhs->parent = move;
      return move;
    }

    /**
     * \brief Creates a move expression of the right hand side tensor
     * expression rhs in the left hand tensor expression lhs.
     * Not indended for direct invocation. Use Move::create instead
     **/
    Move(
      const PTR(IndexedTensor<F>) &lhs_,
      const PTR(Expression<F>) &rhs_,
      const F beta_,
      const typename Expression<F>::ProtectedToken &
    ) {
      static_assert(
        cc4s::StaticAssert<F>::FALSE,
        "Only tensors or contractions may be used as the right hand side of an assignment."
      );
    }
    /**
     * \brief Creates a move expression where the right hand side is
     * single IndexedTensor expression. It will be wrapped in a Contraction
     * Expression so that all moves have contractions as their right hand side.
     * Not indended for direct invocation. Use Move::create instead
     **/
    Move(
      const PTR(IndexedTensor<F>) &lhs_,
      const PTR(IndexedTensor<F>) &rhs_,
      const F beta_,
      const typename Expression<F>::ProtectedToken &
    ): lhs(lhs_), rhs(Contraction<F>::create(1, rhs_)), beta(beta_) {
    }
    /**
     * \brief Creates a move expression where the right hand side is
     * a Contraction expression.
     * Not indended for direct invocation. Use Move::create instead
     **/
    Move(
      const PTR(IndexedTensor<F>) &lhs_,
      const PTR(Contraction<F>) &rhs_,
      const F beta_,
      const typename Expression<F>::ProtectedToken &
    ): lhs(lhs_), rhs(rhs_), beta(beta_) {
    }
    virtual ~Move() {
    }

    virtual void countIndices(IndexCounts &indexCounts) const {
      lhs->countIndices(indexCounts);
      rhs->countIndices(indexCounts);
    }

    PTR(IndexedTensor<F>) lhs;
    PTR(Contraction<F>) rhs;
    F beta;
  };

  /**
   * \brief Creates an updating move expression.
   **/
  template <typename Lhs, typename Rhs>
  inline PTR(Move<typename Lhs::FieldType>) operator +=(
    const PTR(Lhs) &lhs, const PTR(Rhs) &rhs
  ) {
    return Move<typename Lhs::FieldType>::create(lhs, rhs, 1);
  }
  /**
   * \brief Creates an updating move expression.
   **/
  template <typename Lhs, typename Rhs>
  inline PTR(Move<typename Lhs::FieldType>) operator -=(
    const PTR(Lhs) &lhs, const PTR(Rhs) &rhs
  ) {
    return Move<typename Lhs::FieldType>::create(
      lhs, Contraction<typename Rhs::FieldType>::create(-1, rhs),
      1
    );
  }
}

#endif

