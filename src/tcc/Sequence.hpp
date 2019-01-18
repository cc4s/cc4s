/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_SEQUENCE_DEFINED
#define TCC_SEQUENCE_DEFINED

#include <tcc/Expression.hpp>
#include <tcc/OperationSequence.hpp>

#include <util/SharedPointer.hpp>
#include <util/StaticAssert.hpp>

#include <vector>

namespace tcc {
  template <typename F, typename TE>
  class Move;

  template <typename TE>
  class Sequence: public Expression<TE> {
  public:
    /**
     * \brief Creates a sequence expression of the two given tensor
     * expressions A and B.
     **/
    template <typename LHS, typename RHS>
    static PTR(Sequence<typename LHS::TensorEngine>) create(
      const PTR(LHS) &A, const PTR(RHS) &B
    ) {
      static_assert(
        cc4s::TypeRelations<
          typename LHS::TensorEngine, typename RHS::TensorEngine
        >::EQUALS,
        "All expressions within a sequence must have the same tensor engine."
      );
      return NEW(Sequence<typename LHS::TensorEngine>,
        A, B,
        typename Expression<typename LHS::TensorEngine>::ProtectedToken()
      );
    }

    /**
     * \brief Constructor for the empty sequence.
     * Not intended for direct invocation, create sequences
     * using the static create method.
     **/
    // TODO: use protected token
    Sequence() { }

    /**
     * \brief Flattening constructor given two sequences.
     * Not intended for direct invocation, create sequences
     * using the static create method.
     **/
    Sequence(
      const PTR(Sequence<TE>) &lhs,
      const PTR(Sequence<TE>) &rhs,
      const typename Expression<TE>::ProtectedToken &
    ): moves(lhs->moves) {
      moves.insert(moves.end(), rhs->moves.begin(), rhs->moves.end());
    }
    /**
     * \brief Flattening constructor given a sequence on the left hand
     * side and another expression on the right hand side.
     * Not intended for direct invocation, create sequences
     * using the static create method.
     **/
    template <typename F>
    Sequence(
      const PTR(Sequence<TE>) &lhs,
      const PTR(ESC(Move<F,TE>)) &rhs,
      const typename Expression<TE>::ProtectedToken &
    ): moves(lhs->moves) {
      moves.push_back(rhs);
    }
    /**
     * \brief Flattening constructor given a sequence on the right hand
     * side and another expression on the left hand side.
     * Not intended for direct invocation, create sequences
     * using the static create method.
     **/
    template <typename F>
    Sequence(
      const PTR(ESC(Move<F,TE>)) &lhs,
      const PTR(Sequence<TE>) &rhs,
      const typename Expression<TE>::ProtectedToken &
    ) {
      moves.push_back(lhs);
      moves.insert(moves.end(), rhs->moves.begin(), rhs->moves.end());
    }
    /**
     * \brief Constructor given two moves.
     * Not intended for direct invocation, create sequences
     * using the static create method.
     **/
    template <typename F, typename G>
    Sequence(
      const PTR(ESC(Move<F,TE>)) &lhs,
      const PTR(ESC(Move<G,TE>)) &rhs,
      const typename Expression<TE>::ProtectedToken &
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
    template <typename LHS, typename RHS>
    Sequence(
      const PTR(LHS) &lhs,
      const PTR(RHS) &rhs,
      const typename Expression<TE>::ProtectedToken &
    ) {
      static_assert(
        cc4s::StaticAssert<TE>::FALSE,
        "Only sequences of moves are possible."
      );
    }

    virtual ~Sequence() {
    }

    virtual PTR(Operation<TE>) compile(IndexCounts &) {
      std::vector<PTR(Operation<TE>)> operations(moves.size());
      for (size_t i(0); i < moves.size(); ++i) {
        IndexCounts indexCounts;
        operations[i] = moves[i]->compile(indexCounts);
      }
      return OperationSequence<TE>::create(operations);
    }

    virtual void countIndices(IndexCounts &) {
      // the indidex of each subexpression are independet of each other
      // so nothing will be counted.
      // counting will be done on the level of moves and contractions
    }

  protected:
    std::vector<PTR(Expression<TE>)> moves;
  };

  /**
   * \brief Creates a sequence expression of the two given tensor
   * expressions A and B using the comma operator.
   **/
  template <typename LHS, typename RHS>
  inline PTR(Sequence<typename RHS::TensorEngine>) operator ,(
    const PTR(LHS) &A, const PTR(RHS) &B
  ) {
    return Sequence<typename RHS::TensorEngine>::create(A, B);
  }
}

#endif

