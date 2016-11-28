/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_SUM_OPERATION_DEFINED
#define TENSOR_SUM_OPERATION_DEFINED

namespace cc4s {
  template <typename F>
  class TensorSumOperation: public TensorOperation<F> {
  public:
    TensorSumOperation() {
    }
    virtual ~TensorSumOperation() {
    }
    virtual void execute() {
    }
  protected:
  };
}

#endif

