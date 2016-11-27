/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_TENSOR_ASSIGNMENT_DEFINED
#define DRY_TENSOR_ASSIGNMENT_DEFINED

#include <util/DryTensorExpression.hpp>
#include <util/IndexedDryTensor.hpp>

namespace cc4s {
  template <typename F>
  class DryTensorAssignment: public cc4s::DryTensorExpression<F> {
  public:
    DryTensorAssignment(
      IndexedDryTensor<F> *lhs_, DryTensorExpression<F> *rhs_
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
  protected:
    IndexedDryTensor<F> *lhs;
    DryTensorExpression<F> *rhs;
  };
}

#endif

