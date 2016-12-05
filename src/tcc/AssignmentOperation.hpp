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
    AssignmentOperation(
      const std::shared_ptr<FetchOperation<F>> &lhs_,
      const std::shared_ptr<Operation<F>> &rhs_
    ):
      Operation<F>(rhs_->costs),
      lhs(lhs_), rhs(rhs_)
    {
    }
    virtual ~AssignmentOperation() {
    }
    virtual void execute() {
      // TODO: access actual tensors within Tensors for execution
      // lhs->tensor->t->sum(lhs->indices, ... , *rhs->tensor->t, rhs->indices);
      rhs->execute();
      lhs->execute();
      LOG(1, "TCC") << "executing " <<
        lhs->getResult()->get_name() << "[" << lhs->getResultIndices() <<
        "] = " <<
        rhs->getResult()->get_name() << "[" << rhs->getResultIndices() <<
        "]" << std::endl;
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

