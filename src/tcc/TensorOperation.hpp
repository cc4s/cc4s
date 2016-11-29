/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_OPERATION_DEFINED
#define TENSOR_OPERATION_DEFINED

#include <tcc/DryTensor.hpp>
#include <tcc/Costs.hpp>
#include <tcc/TensorExpression.hpp>
#include <util/Log.hpp>
#include <algorithm>

namespace cc4s {
  template <typename F>
  class TensorOperation {
  public:
    TensorOperation(Costs const &costs_): costs(costs_) {
    }
    virtual ~TensorOperation() {
    }

    virtual void execute() = 0;

    virtual DryTensor<F> *getResult() = 0;
    virtual std::string const &getResultIndices() = 0;

    /**
     * \brief Costs to evaluate this operation in time and memory
     **/
    Costs costs;

  protected:
    /**
     * \brief Delete suboperations without deleting fetch operations at the
     * leaves. This is used
     * by the contraction compiler when different orders of contrations
     * need to be tried without deleting the tensor fetch operations each time.
     **/
    virtual void clearLeavingFetches() {
    }
  };

  template <typename F>
  TensorOperation<F> *compile(TensorExpression<F> &expression) {
    return expression.compile("");
  }
}

#endif

