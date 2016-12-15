/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_MOVE_OPERATION_DEFINED
#define TCC_MOVE_OPERATION_DEFINED

#include <tcc/NonVoidOperation.hpp>

#include <memory>

namespace tcc {
  template <typename F>
  class MoveOperation: public NonVoidOperation<F> {
  public:
    /**
     * \brief Creates a move operation moving the results of
     * the right hand side operation into given tensor of the
     * left hand side.
     * Not intended for direct invocation. Use Tcc::compile(expression) to
     * generate operations.
     **/
    MoveOperation(
      const std::shared_ptr<Operation<F>> &rhs_,
      const std::shared_ptr<Tensor<F>> &result_,
      const char *resultIndices_,
      Costs moveCosts,
      const typename Operation<F>::ProtectedToken &
    ):
      NonVoidOperation<F>(
        result_, resultIndices_,
        rhs_->costs,
        typename Operation<F>::ProtectedToken()
      ),
      rhs(rhs_)
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
    static std::shared_ptr<MoveOperation<F>>  create(
      const std::shared_ptr<Operation<F>> &rhs_,
      const std::shared_ptr<Tensor<F>> &result_,
      const char *resultIndices_,
      const Costs &moveCosts
    ) {
      return std::make_shared<MoveOperation<F>>(
        rhs_,
        result_, resultIndices_, moveCosts,
        typename Operation<F>::ProtectedToken()
      );
    }

    std::shared_ptr<Operation<F>> rhs;

    friend class Tcc<F>;
  };
}

#endif

