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
     * \brief Creates a move expression of the right hand side tensor
     * expression rhs in the left hand tensor expression lhs.
     * Not indended for direct invocation. Use Move::create or
     * the operator <<= instead.
     **/
    Move(
      const std::shared_ptr<IndexedTensor<F>> &lhs_,
      const std::shared_ptr<Expression<F>> &rhs_,
      const typename Expression<F>::ProtectedToken &
    ) {
      static_assert(
        cc4s::StaticAssert<F>::False,
        "Only tensors or contractions may be used as the right hand side of an assignment."
      );
    }
    Move(
      const std::shared_ptr<IndexedTensor<F>> &lhs_,
      const std::shared_ptr<IndexedTensor<F>> &rhs_,
      const typename Expression<F>::ProtectedToken &
    ): lhs(lhs_), rhs(rhs_) {
    }
    Move(
      const std::shared_ptr<IndexedTensor<F>> &lhs_,
      const std::shared_ptr<Contraction<F>> &rhs_,
      const typename Expression<F>::ProtectedToken &
    ): lhs(lhs_), rhs(rhs_) {
    }
    virtual ~Move() {
    }

    virtual std::shared_ptr<Operation<F>> compile(const std::string &) {
      return std::make_shared<MoveOperation<F>>(
        std::make_shared<FetchOperation<F>>(
          lhs,
          typename Operation<F>::ProtectedToken()
        ),
        rhs->compile(lhs->indices),
        typename Operation<F>::ProtectedToken()
      );
    }

    /**
     * \brief Moves the given right hand side expression into the given
     * indexed tensor, returning the indexed tensor as result expression
     * for possible further operations.
     **/
    template <typename Rhs>
    static inline std::shared_ptr<Move<typename Rhs::FieldType>> create(
      const std::shared_ptr<IndexedTensor<typename Rhs::FieldType>> &lhs,
      const std::shared_ptr<Rhs> &rhs
    ) {
      return std::make_shared<Move<typename Rhs::FieldType>>(
        lhs, rhs,
        typename Expression<typename Rhs::FieldType>::ProtectedToken()
      );
    }

    std::shared_ptr<IndexedTensor<F>> lhs;
    std::shared_ptr<Expression<F>> rhs;
  };
}

#endif

