/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_SEQUENCE_DEFINED
#define TCC_SEQUENCE_DEFINED

#include <tcc/Expression.hpp>

#include <util/SharedPointer.hpp>
#include <util/StaticAssert.hpp>

#include <vector>

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
    static PTR(Sequence<typename Lhs::FieldType>) create(
      const PTR(Lhs) &A, const PTR(Rhs) &B
    ) {
      static_assert(
        cc4s::TypeRelations<
          typename Lhs::FieldType, typename Rhs::FieldType
        >::EQUALS,
        "Currently, only tensors of the same type can be in a sequence."
      );
      auto sequence(
        NEW(Sequence<typename Lhs::FieldType>,
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
      const PTR(Sequence<F>) &lhs,
      const PTR(Sequence<F>) &rhs,
      const typename Expression<F>::ProtectedToken &
    ): moves(lhs->moves) {
      moves.insert(moves.end(), rhs->moves.begin(), rhs->moves.end());
    }
    /**
     * \brief Flattening constructor given a sequence on the left hand
     * side and another expression on the right hand side.
     * Not intended for direct invocation, create sequences
     * using the static create method.
     **/
    Sequence(
      const PTR(Sequence<F>) &lhs,
      const PTR(Move<F>) &rhs,
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
      const PTR(Move<F>) &lhs,
      const PTR(Sequence<F>) &rhs,
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
      const PTR(Move<F>) &lhs,
      const PTR(Move<F>) &rhs,
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
      const PTR(Expression<F>) &lhs,
      const PTR(Expression<F>) &rhs,
      const typename Expression<F>::ProtectedToken &
    ) {
      static_assert(
        cc4s::StaticAssert<F>::FALSE,
        "Only sequences of moves are supported."
      );
    }

    virtual ~Sequence() {
    }

    std::vector<PTR(Move<F>)> moves;
  };

  /**
   * \brief Creates a sequence expression of the two given tensor
   * expressions A and B using the multiplication operator *.
   **/
  template <typename Lhs, typename Rhs>
  inline PTR(Sequence<typename Lhs::FieldType>) operator ,(
    const PTR(Lhs) &A, const PTR(Rhs) &B
  ) {
    return Sequence<typename Lhs::FieldType>::create(A, B);
  }
}

#endif

