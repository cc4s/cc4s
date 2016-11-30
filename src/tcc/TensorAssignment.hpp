/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_ASSIGNMENT_DEFINED
#define TENSOR_ASSIGNMENT_DEFINED

#include <tcc/TensorExpression.hpp>
#include <tcc/IndexedTensor.hpp>
#include <tcc/TensorOperation.hpp>
#include <tcc/TensorAssignmentOperation.hpp>
#include <util/StaticAssert.hpp>
#include <util/Exception.hpp>

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

    virtual void log() const {
      rhs->log();
      lhs->log();
      LOG(0,"TCC") << "assignment" << std::endl;
    };

    virtual TensorOperation<F> *compile(std::string const &) {
      return new TensorAssignmentOperation<F>(
        new TensorFetchOperation<F>(lhs),
        rhs->compile(lhs->indices)
      );
    }

    IndexedTensor<F> *lhs;
    TensorExpression<F> *rhs;
  };
}

#endif

