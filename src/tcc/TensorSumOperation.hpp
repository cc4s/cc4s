/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_SUM_OPERATION_DEFINED
#define TENSOR_SUM_OPERATION_DEFINED

#include <tcc/IndexedTensor.hpp>

namespace cc4s {
  template <typename F>
  class TensorSumOperation: public TensorOperation<F> {
  public:
    TensorSumOperation(
      IndexedTensor<F> *lhs_, IndexedTensor<F> *rhs_
    ): lhs(lhs_), rhs(rhs_) {
    }
    virtual ~TensorSumOperation() {
    }
    virtual void execute() {
      // TODO: access actual tensors within DryTensors for execution
      // lhs->tensor->t->sum(lhs->indices, ... , *rhs->tensor->t, rhs->indices);
    }

  protected:
    IndexedTensor<F> *lhs, *rhs;
  };
}

#endif

