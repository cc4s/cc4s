/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_TENSOR_ASSIGNMENT_DEFINED
#define DRY_TENSOR_ASSIGNMENT_DEFINED

#include <tcc/DryTensorExpression.hpp>
#include <tcc/IndexedDryTensor.hpp>

namespace cc4s {
  template <typename F>
  class DryTensorAssignment: public cc4s::DryTensorExpression<F> {
  public:
    DryTensorAssignment(
      IndexedDryTensor<F> *lhs_, DryTensorExpression<F> *rhs_
    ) {
      static_assert(
        false, "Only tensors or contractions may be used as the right hand side of an assignment."
      );
    }
    DryTensorAssignment(
      IndexedDryTensor<F> *lhs_, IndexedDryTensor<F> *rhs_
    ): lhs(lhs_), rhs(rhs_) {
    }
    DryTensorAssignment(
      IndexedDryTensor<F> *lhs_, DryTensorContraction<F> *rhs_
    ): lhs(lhs_), rhs(rhs_) {
    }
    virtual ~DryTensorAssignment() {
      delete lhs, delete rhs;
    }
    virtual void log() const {
      rhs->log();
      lhs->log();
      LOG(0,"TCC") << "assignment" << std::endl;
    };

    IndexedDryTensor<F> *lhs;
    DryTensorExpression<F> *rhs;
  };
}

#endif

