/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_OPERATION_DEFINED
#define TENSOR_OPERATION_DEFINED

#include <tcc/DryTensor.hpp>
#include <tcc/Costs.hpp>
#include <tcc/TensorExpression.hpp>
#include <util/Log.hpp>
#include <algorithm>

#include <memory>
using std::shared_ptr;

namespace cc4s {
  template <typename F>
  class TensorContractionOperation;

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
  };

  template <typename F>
  shared_ptr<TensorOperation<F>> compile(TensorExpression<F> &expression) {
    shared_ptr<TensorOperation<F>> result(compile(&expression));
    // TODO: use shared pointers in entire tcc
    delete &expression;
    return result;
  }

  template <typename F>
  shared_ptr<TensorOperation<F>> compile(TensorExpression<F> *expression) {
    return expression->compile("");
  }
}

#endif

