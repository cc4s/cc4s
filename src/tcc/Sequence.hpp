/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_SEQUENCE_DEFINED
#define TCC_SEQUENCE_DEFINED

#include <tcc/Expression.hpp>
#include <util/StaticAssert.hpp>

#include <vector>
#include <memory>

namespace tcc {
  template <typename F>
  class Move;

  template <typename F>
  class Sequence: public Expression<F> {
  public:
    /**
     * \brief Creates a sequence expression of the two given tensor
     * expressions A and B.
     **/
    template <typename Lhs, typename Rhs>
    static std::shared_ptr<Sequence<typename Lhs::FieldType>> create(
      const std::shared_ptr<Lhs> &A, const std::shared_ptr<Rhs> &B
    ) {
      static_assert(
        cc4s::TypeRelations<
          typename Lhs::FieldType, typename Rhs::FieldType
        >::EQUALS,
        "Currently, only tensors of the same type can be in a sequence."
      );
      auto sequence(
        std::make_shared<Sequence<typename Lhs::FieldType>>(
          A, B,
          typename Expression<typename Lhs::FieldType>::ProtectedToken()
        )
      );
      A->parent = sequence;
      B->parent = sequence;
      return sequence;
    }

    /**
     * \brief Flattening constructor given two sequences.
     * Not intended for direct invocation, create sequences
     * using the static create method.
     **/
    Sequence(
      const std::shared_ptr<Sequence<F>> &lhs,
      const std::shared_ptr<Sequence<F>> &rhs,
      const typename Expression<F>::ProtectedToken &
    ): alpha(lhs->alpha * rhs->alpha), factors(lhs->factors) {
      factors.insert(factors.end(), rhs->factors.begin(), rhs->factors.end());
    }
    /**
     * \brief Flattening constructor given a sequence on the left hand
     * side and another expression on the right hand side.
     * Not intended for direct invocation, create sequences
     * using the static create method.
     **/
    Sequence(
      const std::shared_ptr<Sequence<F>> &lhs,
      const std::shared_ptr<Move<F>> &rhs,
      const typename Expression<F>::ProtectedToken &
    ): moves(lhs->moves) {
      moves.push_back(rhs);
    }
    /**
     * \brief Flattening constructor given a sequence on the right hand
     * side and another expression on the left hand side.
     * Not intended for direct invocation, create sequences
     * using the static create method.
     **/
    Sequence(
      const std::shared_ptr<Move<F>> &lhs,
      const std::shared_ptr<Sequence<F>> &rhs,
      const typename Expression<F>::ProtectedToken &
    ): moves(rhs->moves) {
      moves.push_back(lhs);
    }
    /**
     * \brief Constructor given two moves.
     * Not intended for direct invocation, create sequences
     * using the static create method.
     **/
    Sequence(
      const std::shared_ptr<Move<F>> &lhs,
      const std::shared_ptr<Move<F>> &rhs,
      const typename Expression<F>::ProtectedToken &
    ) {
      moves.push_back(lhs);
      moves.push_back(rhs);
    }
    /**
     * \brief Constructor given two general expressions.
     * This is currently not supported.
     * Not intended for direct invocation, create sequences
     * using the static create method.
     **/
    Sequence(
      const std::shared_ptr<Expression<F>> &lhs,
      const std::shared_ptr<Expression<F>> &rhs,
      const typename Expression<F>::ProtectedToken &
    ) {
      static_assert(
        cc4s::StaticAssert<F>::FALSE,
        "Only sequences of moves are supported."
      );
    }

    virtual ~Sequence() {
    }

    std::vector<std::shared_ptr<Move<F>>> moves;
  };

  /**
   * \brief Creates a sequence expression of the two given tensor
   * expressions A and B using the multiplication operator *.
   **/
  template <typename Lhs, typename Rhs>
  inline std::shared_ptr<Sequence<typename Lhs::FieldType>> operator ,(
    const std::shared_ptr<Lhs> &A, const std::shared_ptr<Rhs> &B
  ) {
    return Sequence<typename Lhs::FieldType>::create(A, B);
  }
}

#endif

