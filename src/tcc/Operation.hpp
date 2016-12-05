/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_OPERATION_DEFINED
#define TCC_OPERATION_DEFINED

#include <tcc/Tensor.hpp>
#include <tcc/Costs.hpp>
#include <tcc/Expression.hpp>
#include <util/Log.hpp>
#include <algorithm>

#include <memory>

namespace tcc {
  template <typename F>
  class ContractionOperation;

  template <typename F>
  class Operation {
  public:
    Operation(const Costs &costs_): costs(costs_) {
    }
    virtual ~Operation() {
    }

    virtual void execute() = 0;

    virtual std::shared_ptr<Tensor<F>> getResult() = 0;
    virtual std::string const &getResultIndices() = 0;

    /**
     * \brief Costs to evaluate this operation in time and memory
     **/
    Costs costs;
  };

  template <typename F>
  std::shared_ptr<Operation<F>> compile(
    const std::shared_ptr<Expression<F>> &expression
  ) {
    return compile(expression, "");
  }
}

#endif

