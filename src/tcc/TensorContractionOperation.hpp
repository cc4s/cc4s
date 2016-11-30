/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_CONTRACTION_OPERATION_DEFINED
#define TENSOR_CONTRACTION_OPERATION_DEFINED

namespace cc4s {
  template <typename F>
  class TensorContractionOperation: public TensorOperation<F> {
  public:
    TensorContractionOperation() {
    }
    virtual ~TensorContractionOperation() {
    }
    virtual void execute() {
    }
  protected:
  };
}

#endif

