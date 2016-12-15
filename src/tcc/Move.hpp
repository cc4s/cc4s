/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_ASSIGNMENT_DEFINED
#define TCC_ASSIGNMENT_DEFINED

#include <tcc/Expression.hpp>
#include <tcc/MoveOperation.hpp>
#include <tcc/FetchOperation.hpp>
#include <util/StaticAssert.hpp>
#include <util/Exception.hpp>

#include <memory>

namespace tcc {
  template <typename F>
  class Move: public Expression<F> {
  public:
    /**
     * \brief Moves the given right hand side expression into the given
     * indexed tensor, returning the indexed tensor as result expression
     * for possible further operations.
     **/
    template <typename Lhs, typename Rhs, typename S>
    static inline std::shared_ptr<Move<typename Lhs::FieldType>> create(
      const std::shared_ptr<Lhs> &lhs,
      const std::shared_ptr<Rhs> &rhs,
      const S beta
    ) {
      static_assert(
        cc4s::TypeRelations<S, typename Lhs::FieldType>::CASTABLE_TO,
        "The type of the scalar must be convertible to the tensor type."
      );
      auto move(
        std::make_shared<Move<typename Lhs::FieldType>>(
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
      const std::shared_ptr<IndexedTensor<F>> &lhs_,
      const std::shared_ptr<Expression<F>> &rhs_,
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
      const std::shared_ptr<IndexedTensor<F>> &lhs_,
      const std::shared_ptr<IndexedTensor<F>> &rhs_,
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
      const std::shared_ptr<IndexedTensor<F>> &lhs_,
      const std::shared_ptr<Contraction<F>> &rhs_,
      const F beta_,
      const typename Expression<F>::ProtectedToken &
    ): lhs(lhs_), rhs(rhs_), beta(beta_) {
    }
    virtual ~Move() {
    }

    std::shared_ptr<IndexedTensor<F>> lhs;
    std::shared_ptr<Contraction<F>> rhs;
    F beta;
  };

  /**
   * \brief Creates an updating move expression.
   **/
  template <typename Lhs, typename Rhs>
  inline std::shared_ptr<Move<typename Lhs::FieldType>> operator +=(
    const std::shared_ptr<Lhs> &lhs, const std::shared_ptr<Rhs> &rhs
  ) {
    return Move<typename Lhs::FieldType>::create(lhs, rhs, 1);
  }
  /**
   * \brief Creates an updating move expression.
   **/
  template <typename Lhs, typename Rhs>
  inline std::shared_ptr<Move<typename Lhs::FieldType>> operator -=(
    const std::shared_ptr<Lhs> &lhs, const std::shared_ptr<Rhs> &rhs
  ) {
    return Move<typename Lhs::FieldType>::create(
      lhs, Contraction<typename Rhs::FieldType>::create(-1, rhs),
      1
    );
  }
}

#endif

