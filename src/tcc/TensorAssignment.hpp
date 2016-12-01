/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_ASSIGNMENT_DEFINED
#define TENSOR_ASSIGNMENT_DEFINED

#include <tcc/TensorExpression.hpp>
#include <tcc/IndexedTensor.hpp>
#include <tcc/TensorOperation.hpp>
#include <tcc/TensorAssignmentOperation.hpp>
#include <util/StaticAssert.hpp>
#include <util/Exception.hpp>

#include <memory>
using std::shared_ptr;
using std::make_shared;

namespace cc4s {
  template <typename F>
  class TensorAssignment: public TensorExpression<F> {
  public:
    TensorAssignment(
      IndexedTensor<F> *lhs_, TensorExpression<F> *rhs_
    ) {
      static_assert(
        StaticAssert<F>::False,
        "Only tensors or contractions may be used as the right hand side of an assignment."
      );
    }
    TensorAssignment(
      IndexedTensor<F> *lhs_, IndexedTensor<F> *rhs_
    ): lhs(lhs_), rhs(rhs_) {
    }
    TensorAssignment(
      IndexedTensor<F> *lhs_, TensorContraction<F> *rhs_
    ): lhs(lhs_), rhs(rhs_) {
    }
    virtual ~TensorAssignment() {
      delete lhs, delete rhs;
    }

    virtual shared_ptr<TensorOperation<F>> compile(std::string const &) {
      return make_shared<TensorAssignmentOperation<F>>(
        make_shared<TensorFetchOperation<F>>(lhs),
        rhs->compile(lhs->indices)
      );
    }

    IndexedTensor<F> *lhs;
    TensorExpression<F> *rhs;
  };
}

#endif

