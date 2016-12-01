/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_ASSIGNMENT_OPERATION_DEFINED
#define TENSOR_ASSIGNMENT_OPERATION_DEFINED

#include <tcc/TensorOperation.hpp>
#include <tcc/TensorFetchOperation.hpp>

namespace cc4s {
  template <typename F>
  class TensorAssignmentOperation: public TensorOperation<F> {
  public:
    TensorAssignmentOperation(
      TensorFetchOperation<F> *lhs_, TensorOperation<F> *rhs_
    ):
      TensorOperation<F>(rhs_->costs),
      lhs(lhs_), rhs(rhs_)
    {
    }
    virtual ~TensorAssignmentOperation() {
      // the suboperations are dependent entities
      if (lhs) delete lhs;
      if (rhs) delete rhs;
    }
    virtual void execute() {
      // TODO: access actual tensors within DryTensors for execution
      // lhs->tensor->t->sum(lhs->indices, ... , *rhs->tensor->t, rhs->indices);
      rhs->execute();
      lhs->execute();
      LOG(0, "TCC") << "executing " <<
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
    TensorFetchOperation<F> *lhs;
    TensorOperation<F> *rhs;
  };
}

#endif

