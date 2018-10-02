/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_OPERATION_DEFINED
#define TCC_OPERATION_DEFINED

#include <tcc/Costs.hpp>

#include <util/SharedPointer.hpp>

namespace tcc {
  template <typename F>
  class Tensor;

  template <typename F>
  class Operation {
  public:
    virtual ~Operation() {
    }

    virtual void execute() = 0;

    virtual PTR(Tensor<F>) getResult() = 0;
    virtual const std::string &getResultIndices() = 0;

    /**
     * \brief Costs to evaluate this operation in time and memory
     **/
    Costs costs;

  protected:
    Operation(const Costs &costs_): costs(costs_) {
    }
    class ProtectedToken {
    };

    friend class Tcc<F>;
  };
}

#endif

