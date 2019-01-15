/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_MOVE_DEFINED
#define TCC_MOVE_DEFINED

#include <tcc/Expression.hpp>
#include <tcc/Contraction.hpp>
#include <tcc/MoveOperation.hpp>

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
/*
      lhs->parent = move;
      rhs->parent = move;
*/
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

    // each move has its private index namespace so disregard the outer
    // indexCounts
    virtual PTR(Operation<F>) compile(IndexCounts &) {
      LOG(2, "TCC") << "compiling move..." << std::endl;
      // create a new namespace of indices
      IndexCounts indexCounts;
      // and determine how often each index is used
      countIndices(indexCounts);

      // TODO: currently only move(lhs,contraction(factors...)) are done
      indexCounts.triedPossibilitiesCount = 0;
      auto operation(
        // compile right-hand-side in the namespace of this move
        DYNAMIC_PTR_CAST(TensorResultOperation<F>, rhs->compile(indexCounts))
      );
      // write operation results directly to lhs tensor instead of intermediate
      operation->result = lhs->tensor;
      operation->resultIndices = lhs->indices;
      // enter the beta factor of this move
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
    PTR(IndexedTensor<F>) lhs;
    PTR(Contraction<F>) rhs;
    F beta;

    friend class Tcc<F>;
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

  /**
   * \brief Moves the given right hand side expression into the given
   * indexed tensor, returning the indexed tensor as result expression
   * for possible further operations.
   * Note that the operator = cannot be used since all expressions are
   * represented by shared pointers.
   **/
  template <typename Rhs>
  inline PTR(Move<typename Rhs::FieldType>) operator <<=(
    const PTR(IndexedTensor<typename Rhs::FieldType>) &lhs,
    const PTR(Rhs) &rhs
  ) {
    return Move<typename Rhs::FieldType>::create(lhs, rhs, 0);
  }

  template <typename Lhs, typename Rhs>
  inline PTR(Move<typename Lhs::FieldType>) operator <<=(
    const PTR(Lhs) &, const PTR(Rhs) &rhs
  ) {
    static_assert(
      cc4s::TypeRelations<
        typename Lhs::FieldType, typename Rhs::FieldType
      >::EQUALS,
      "Move operations requires tensors of same type."
    );
    static_assert(
      cc4s::StaticAssert<Lhs>::FALSE,
      "Only indexed tensors may be used as the left hand side of a move operation."
    );
    return PTR(Move<typename Lhs::FieldType>)();
  }

  template <typename F>
  inline std::ostream &operator <<(
    std::ostream &stream, const IndexedTensor<F> &t
  ) {
    return stream << t.tensor->getName() << "[" << t.indices << "]";
  }
}

#endif

