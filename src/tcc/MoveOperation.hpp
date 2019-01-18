/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_MOVE_OPERATION_DEFINED
#define TCC_MOVE_OPERATION_DEFINED

#include <tcc/TensorResultOperation.hpp>

#include <util/SharedPointer.hpp>

namespace tcc {
  template <typename F, typename TE> class Contraction;

  template <typename F, typename TE>
  class MoveOperation: public TensorResultOperation<F,TE> {
  public:
    /**
     * \brief Creates a move operation moving the results of
     * the right hand side operation into given tensor of the
     * left hand side after applying the function f.
     * The function f defaults to the identity operation.
     * Not intended for direct invocation. Use Tcc::compile(expression) to
     * generate operations.
     **/
    MoveOperation(
      const PTR(ESC(TensorResultOperation<F,TE>)) &rhs_,
      const PTR(ESC(Tensor<F,TE>)) &result_,
      const char *resultIndices_,
      Costs moveCosts,
      const typename Operation<TE>::ProtectedToken &
    ):
      TensorResultOperation<F,TE>(
        result_, resultIndices_,
        rhs_->costs,
        typename Operation<TE>::ProtectedToken()
      ),
      rhs(rhs_)
    {
      this->costs += moveCosts;
    }

    virtual ~MoveOperation() {
    }

    virtual void execute() {
      rhs->execute();
      this->getResult()->getMachineTensor()->move(
        this->alpha,
        rhs->getResult()->getMachineTensor(), rhs->getResultIndices(),
        this->beta,
        this->resultIndices
      );
    }

  protected:
    static PTR(ESC(MoveOperation<F,TE>)) create(
      const PTR(ESC(TensorResultOperation<F,TE>)) &rhs_,
      const PTR(ESC(Tensor<F,TE>)) &result_,
      const char *resultIndices_,
      const Costs &moveCosts
    ) {
      return NEW(ESC(MoveOperation<F,TE>),
        rhs_,
        result_, resultIndices_, moveCosts,
        typename Operation<TE>::ProtectedToken()
      );
    }

    PTR(ESC(TensorResultOperation<F,TE>)) rhs;
    F alpha, beta;

    friend class Contraction<F,TE>;
  };
}

#endif

