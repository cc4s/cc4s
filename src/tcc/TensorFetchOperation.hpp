/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_FETCH_OPERATION_DEFINED
#define TENSOR_FETCH_OPERATION_DEFINED

#include <tcc/DryTensor.hpp>
#include <tcc/IndexedTensor.hpp>

namespace cc4s {
  template <typename F>
  class TensorFetchOperation: public TensorOperation<F> {
  public:
    TensorFetchOperation(
      IndexedTensor<F> *t_
    ):
      TensorOperation<F>(Costs(t_->tensor->getElementsCount())),
      tensor(t_->tensor),
      indices(t_->indices)
    {
    }
    virtual ~TensorFetchOperation() {
      // the fetched tensor is an idependent entity and not to be deleted
    }

    virtual void execute() {
      // nothing to be done in a fetch
    }

    virtual DryTensor<F> *getResult() {
      return tensor;
    }

    virtual std::string const &getResultIndices() {
      return indices;
    }

  protected:
    DryTensor<F> *tensor;
    std::string indices;
  };
}

#endif

