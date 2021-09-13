#ifndef TCC_SEQUENCE_DEFINED
#define TCC_SEQUENCE_DEFINED

#include <tcc/Expression.hpp>

#include <tcc/SequenceOperation.hpp>
#include <util/SharedPointer.hpp>
#include <util/StaticAssert.hpp>
#include <vector>

namespace cc4s {
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
        TypeRelations<
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
        StaticAssert<TE>::FALSE,
        "Only sequences of moves are possible."
      );
    }

    PTR(Operation<TE>) compile(Scope &outerScope) override {
      std::vector<PTR(Operation<TE>)> operations(moves.size());
      for (size_t i(0); i < moves.size(); ++i) {
        operations[i] = moves[i]->compile(outerScope);
      }
      return SequenceOperation<TE>::create(operations, outerScope);
    }

    PTR(Operation<TE>) compile(
      const std::string &file, const size_t line
    ) override {
      Scope scope(file, line);
      return this->compile(scope);
    }


    virtual void countIndices(Scope &) {
      // the indidex of each subexpression are independet of each other
      // so nothing will be counted.
      // counting will be done on the level of moves and contractions
    }

    virtual operator std::string () const {
      std::stringstream stream;
      stream << "Sequence( ";
      std::string delimiter("");
      for (auto const &move: moves) {
        stream << delimiter << std::string(*move);
        delimiter = ", ";
      }
      stream << " )";
      return stream.str();
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

