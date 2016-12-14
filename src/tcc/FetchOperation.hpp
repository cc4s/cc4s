/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_FETCH_OPERATION_DEFINED
#define TCC_FETCH_OPERATION_DEFINED

#include <tcc/Tensor.hpp>
#include <util/Log.hpp>

#include <memory>

namespace tcc {
  template <typename F>
  class FetchOperation: public Operation<F> {
  public:
    /**
     * \brief Creates a fetch operation of a tensor making it accessible
     * for subsequent move or contraction operations.
     * Not intended for direct invocation. Use compile(expression) to
     * generate operations.
     **/
    FetchOperation(
      const std::shared_ptr<Tensor<F>> &tensor_,
      const std::string &indices_,
      const typename Operation<F>::ProtectedToken &
    ):
      Operation<F>(Costs(tensor_->getElementsCount())),
      tensor(tensor_),
      indices(indices_)
    {
    }
    virtual ~FetchOperation() {
    }

    virtual void execute() {
    }

    virtual std::shared_ptr<Tensor<F>> getResult() {
      return tensor;
    }

    virtual const std::string &getResultIndices() {
      return indices;
    }

  protected:
    std::shared_ptr<Tensor<F>> tensor;
    std::string indices;
  };
}

#endif

