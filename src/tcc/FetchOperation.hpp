/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_FETCH_OPERATION_DEFINED
#define TCC_FETCH_OPERATION_DEFINED

#include <tcc/Operation.hpp>
#include <tcc/Costs.hpp>

#include <memory>

namespace tcc {
  template <typename F> class Tensor;
  template <typename F> class IndexedTensor;

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
      const PTR(IndexedTensor<F>) &t,
      const typename Operation<F>::ProtectedToken &
    ):
      Operation<F>(Costs(t->getTensor()->getElementsCount())),
      tensor(t->getTensor()),
      indices(t->getIndices())
    {
    }
    virtual ~FetchOperation() {
    }

    virtual void execute() {
    }

    virtual PTR(Tensor<F>) getResult() {
      return tensor;
    }

    virtual const std::string &getResultIndices() {
      return indices;
    }

  protected:
    static PTR(FetchOperation<F>) create(
      const PTR(IndexedTensor<F>) &indexedTensor
    ) {
      return NEW(FetchOperation<F>,
        indexedTensor, typename Operation<F>::ProtectedToken()
      );
    }

    PTR(Tensor<F>) tensor;
    std::string indices;

    friend class IndexedTensor<F>;
  };
}

#endif

