/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_ASSIGNMENT_OPERATION_DEFINED
#define TENSOR_ASSIGNMENT_OPERATION_DEFINED

#include <tcc/TensorOperation.hpp>
#include <tcc/TensorFetchOperation.hpp>

#include <memory>
using std::shared_ptr;

namespace cc4s {
  template <typename F>
  class TensorAssignmentOperation: public TensorOperation<F> {
  public:
    TensorAssignmentOperation(
      const shared_ptr<TensorFetchOperation<F>> &lhs_,
      const shared_ptr<TensorOperation<F>> &rhs_
    ):
      TensorOperation<F>(rhs_->costs),
      lhs(lhs_), rhs(rhs_)
    {
    }
    virtual ~TensorAssignmentOperation() {
    }
    virtual void execute() {
      // TODO: access actual tensors within DryTensors for execution
      // lhs->tensor->t->sum(lhs->indices, ... , *rhs->tensor->t, rhs->indices);
      rhs->execute();
      lhs->execute();
      LOG(1, "TCC") << "executing " <<
        lhs->getResult()->get_name() << "[" << lhs->getResultIndices() <<
        "] = " <<
        rhs->getResult()->get_name() << "[" << rhs->getResultIndices() <<
        "]" << std::endl;
    }

    virtual DryTensor<F> *getResult() {
      return lhs->getResult();
    }

    virtual std::string const &getResultIndices() {
      return lhs->getResultIndices();
    }

  protected:
    shared_ptr<TensorFetchOperation<F>> lhs;
    shared_ptr<TensorOperation<F>> rhs;
  };
}

#endif

