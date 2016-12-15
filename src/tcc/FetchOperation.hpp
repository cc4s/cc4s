/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_FETCH_OPERATION_DEFINED
#define TCC_FETCH_OPERATION_DEFINED

#include <tcc/Operation.hpp>
#include <tcc/Tensor.hpp>
#include <tcc/Costs.hpp>

#include <memory>

namespace tcc {
  template <typename F>
  class FetchOperation: public Operation<F> {
  public:
    /**
     * \brief Creates a fetch operation of a tensor making it accessible
     * for subsequent move or contraction operations.
     * Not intended for direct invocation. Use Tcc::compile(expression) to
     * generate operations.
     **/
    FetchOperation(
      const std::shared_ptr<IndexedTensor<F>> &t,
      const typename Operation<F>::ProtectedToken &
    ):
      Operation<F>(Costs(t->tensor->getElementsCount())),
      tensor(t->tensor),
      indices(t->indices)
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
    static std::shared_ptr<FetchOperation<F>> create(
      const std::shared_ptr<IndexedTensor<F>> &indexedTensor
    ) {
      return std::make_shared<FetchOperation<F>>(
        indexedTensor, ProtectedToken()
      );
    }

    std::shared_ptr<Tensor<F>> tensor;
    std::string indices;

    friend class Tcc<F>;
  };
}

#endif

