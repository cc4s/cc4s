/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_MOVE_OPERATION_DEFINED
#define TCC_MOVE_OPERATION_DEFINED

#include <tcc/TensorResultOperation.hpp>

#include <util/SharedPointer.hpp>

namespace tcc {
  template <typename F>
  class MoveOperation: public TensorResultOperation<F> {
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
      const PTR(Operation<F>) &rhs_,
      const PTR(Tensor<F>) &result_,
      const char *resultIndices_,
      Costs moveCosts,
      const typename Operation<F>::ProtectedToken &
    ):
      TensorResultOperation<F>(
        result_, resultIndices_,
        rhs_->costs,
        typename Operation<F>::ProtectedToken()
      ),
      rhs(rhs_),
      f(nullptr)
    {
      this->costs += moveCosts;
    }

    MoveOperation(
      const PTR(Operation<F>) &rhs_,
      const std::function<F(const F)> &f_,
      const PTR(Tensor<F>) &result_,
      const char *resultIndices_,
      Costs moveCosts,
      const typename Operation<F>::ProtectedToken &
    ):
      TensorResultOperation<F>(
        result_, resultIndices_,
        rhs_->costs,
        typename Operation<F>::ProtectedToken()
      ),
      rhs(rhs_),
      f(f_)
    {
      this->costs += moveCosts;
    }

    virtual ~MoveOperation() {
    }

    virtual void execute() {
      rhs->execute();
      this->result->getMachineTensor()->move(
        this->alpha,
        rhs->getResult()->getMachineTensor(), rhs->getResultIndices(),
        this->beta,
        this->resultIndices
      );
    }

  protected:
    static PTR(MoveOperation<F>)  create(
      const PTR(Operation<F>) &rhs_,
      const PTR(Tensor<F>) &result_,
      const char *resultIndices_,
      const Costs &moveCosts
    ) {
      return NEW(MoveOperation<F>,
        rhs_,
        result_, resultIndices_, moveCosts,
        typename Operation<F>::ProtectedToken()
      );
    }

    PTR(Operation<F>) rhs;
    std::function<F(const F)> f;

    friend class Tcc<F>;
  };
}

#endif

