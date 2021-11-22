/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef TCC_MOVE_DEFINED
#define TCC_MOVE_DEFINED

#include <tcc/IndexedTensorExpression.hpp>

#include <tcc/Contraction.hpp>
#include <util/SharedPointer.hpp>
#include <util/StaticAssert.hpp>

namespace cc4s {
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
    Ptr<Move<typename LHS::FieldType, typename LHS::TensorEngine>> create(
      const Ptr<LHS> &lhs,
      const Ptr<RHS> &rhs,
      const S beta
    ) {
      static_assert(
        TypeRelations<S, typename LHS::FieldType>::CASTABLE_TO,
        "The type of the scalar must be convertible to the tensor type."
      );
      return New<Move<typename LHS::FieldType,TE>>(
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
      const Ptr<IndexedTensor<F,TE>> &lhs_,
      const Ptr<Expression<TE>> &rhs_,
      const F beta_,
      const typename Expression<TE>::ProtectedToken &
    ) {
      static_assert(
        StaticAssert<F>::FALSE,
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
      const Ptr<Indexing<F,TE>> &lhs_,
      const Ptr<IndexedTensorExpression<F,TE>> &rhs_,
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
      const Ptr<Indexing<F,TE>> &lhs_,
      const Ptr<Contraction<F,TE>> &rhs_,
      const F beta_,
      const typename Expression<TE>::ProtectedToken &
    ): lhs(lhs_), rhs(rhs_), beta(beta_) {
    }
    virtual ~Move() {
    }

    // each move has its private index namespace so disregard the outer
    // scope
    Ptr<Operation<TE>> compile(Scope &outerScope) override {
      LOG_LOCATION(SourceLocation(outerScope.file, outerScope.line)) <<
        "compiling: " << static_cast<std::string>(*this) << std::endl;

      // create a new namespace of indices
      Scope scope(outerScope.file, outerScope.line);
      // and determine how often each index is used
      countIndices(scope);

      // NOTE: currently only move(lhs,contraction(factors...)) are done
      scope.triedPossibilitiesCount = 0;
      auto operation(
        // compile right-hand-side in the namespace of this move
        dynamicPtrCast<IndexedTensorOperation<F,TE>>(rhs->compile(scope))
      );

      LOG_LOCATION(SourceLocation(outerScope.file, outerScope.line))
        << "possibilites tried=" << scope.triedPossibilitiesCount
        << ", best has " << std::string(operation->costs)
        << ": " << std::string(*operation)
        << std::endl;

      operation->beta = beta;
      // compile writing the result to the left-hand-side
      return lhs->lhsCompile(operation);
    }

    // keep other overloads visible
    using Expression<TE>::compile;

    void countIndices(Scope &scope) override {
      lhs->countIndices(scope);
      rhs->countIndices(scope);
    }

    operator std::string () const override {
      std::stringstream stream;
      stream << "Sum( " <<
        std::string(*lhs) << ", " << std::string(*rhs) << ", beta=" << beta << " )";
      return stream.str();
    }

  protected:
    Ptr<IndexedTensorExpression<F,TE>> lhs;
    Ptr<Contraction<F,TE>> rhs;
    F beta;
  };

  /**
   * \brief Creates an updating move expression.
   **/
  template <typename RHS>
  inline
  Ptr<Move<typename RHS::FieldType,typename RHS::TensorEngine>>
  operator +=(
    const
    Ptr<
      Indexing<typename RHS::FieldType,typename RHS::TensorEngine>
    > &lhs,
    const Ptr<RHS> &rhs
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
  Ptr<Move<typename RHS::FieldType,typename RHS::TensorEngine>>
  operator -=(
    const
    Ptr<
      Indexing<typename RHS::FieldType,typename RHS::TensorEngine>
    > &lhs,
    const Ptr<RHS> &rhs
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
  Ptr<Move<typename RHS::FieldType, typename RHS::TensorEngine>>
  operator <<=(
    const
    Ptr<
      Indexing<typename RHS::FieldType,typename RHS::TensorEngine>
    > &lhs,
    const Ptr<RHS> &rhs
  ) {
    return Move<typename RHS::FieldType, typename RHS::TensorEngine>::create(
      lhs, rhs, typename RHS::FieldType(0)
    );
  }

  template <typename LHS, typename RHS>
  inline
  Ptr<Move<typename RHS::FieldType, typename RHS::TensorEngine>>
  operator <<=(
    const Ptr<LHS> &, const Ptr<RHS> &rhs
  ) {
    static_assert(
      TypeRelations<
        typename LHS::FieldType, typename RHS::FieldType
      >::EQUALS,
      "Move operations requires tensors of same type."
    );
    static_assert(
      StaticAssert<RHS>::FALSE,
      "Only indexed tensors may be used as the left hand side of a move operation."
    );
    return Ptr<Move<typename RHS::FieldType,typename RHS::TensorEngine>>();
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

