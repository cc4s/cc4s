/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_ASSIGNMENT_OPERATION_DEFINED
#define TCC_ASSIGNMENT_OPERATION_DEFINED

#include <tcc/Operation.hpp>
#include <tcc/FetchOperation.hpp>

#include <memory>

namespace tcc {
  template <typename F>
  class AssignmentOperation: public Operation<F> {
  public:
    /**
     * \brief Creates a move operation moving the results of
     * the right hand side operation into the fetched tensor on the
     * left hand side.
     * Not intended for direct invocation. Use compile(expression) to
     * generate operations.
     **/
    AssignmentOperation(
      const std::shared_ptr<FetchOperation<F>> &lhs_,
      const std::shared_ptr<Operation<F>> &rhs_,
      const typename Operation<F>::ProtectedToken &
    ):
      Operation<F>(rhs_->costs),
      lhs(lhs_), rhs(rhs_)
    {
    }

    virtual ~AssignmentOperation() {
    }

    virtual void execute() {
      rhs->execute();
      lhs->execute();
      // TODO: access actual tensors within Tensors for execution
      // lhs->tensor->t->sum(lhs->indices, ... , *rhs->tensor->t, rhs->indices);
    }

    virtual std::shared_ptr<Tensor<F>> getResult() {
      return lhs->getResult();
    }

    virtual std::string const &getResultIndices() {
      return lhs->getResultIndices();
    }

  protected:
    std::shared_ptr<FetchOperation<F>> lhs;
    std::shared_ptr<Operation<F>> rhs;
  };
}

#endif

