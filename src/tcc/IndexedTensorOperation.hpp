/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_FETCH_OPERATION_DEFINED
#define TCC_FETCH_OPERATION_DEFINED

#include <tcc/TensorResultOperation.hpp>
#include <tcc/Costs.hpp>

#include <util/SharedPointer.hpp>

namespace tcc {
  template <typename F, typename TE> class IndexedTensor;

  template <typename F, typename TE>
  class FetchOperation: public TensorResultOperation<F,TE> {
  public:
    /**
     * \brief Creates a fetch operation of a tensor making it accessible
     * for subsequent move or contraction operations.
     * Not intended for direct invocation. expression->compile() to
     * generate operations.
     **/
    FetchOperation(
      const PTR(ESC(Tensor<F,TE>)) &t,
      const std::string &indices,
      const typename Operation<TE>::ProtectedToken &
    ):
      TensorResultOperation<F,TE>(
        t, indices.c_str(),
        Costs(t->getElementsCount()),
        typename Operation<TE>::ProtectedToken()
      )
    {
    }
    virtual ~FetchOperation() {
    }

    virtual void execute() {
    }

  protected:
    static PTR(ESC(FetchOperation<F,TE>)) create(
      const PTR(ESC(Tensor<F,TE>)) &t,
      const std::string &indices
    ) {
      return NEW(ESC(FetchOperation<F,TE>),
        t, indices, typename Operation<TE>::ProtectedToken()
      );
    }

    friend class IndexedTensor<F,TE>;
  };
}

#endif

