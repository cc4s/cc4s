/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_NON_VOID_OPERATION_DEFINED
#define TCC_NON_VOID_OPERATION_DEFINED

#include <tcc/Operation.hpp>
#include <tcc/Costs.hpp>
#include <util/Log.hpp>

#include <string>
#include <memory>

namespace tcc {
  template <typename F>
  class Tensor;

  template <typename F>
  class NonVoidOperation: public Operation<F> {
  public:
    /**
     * \brief Creates a contraction operation contracting the results of
     * the left and the right sub-operations where the result is to be
     * stored in the specified result tensor.
     * Not intended for direct invocation. Use Tcc::compile(expression) to
     * generate operations.
     **/
    NonVoidOperation(
      const std::shared_ptr<Tensor<F>> &result_,
      const char *resultIndices_,
      const Costs &costs_,
      const typename Operation<F>::ProtectedToken &
    ):
      Operation<F>(costs_),
      alpha(F(1)), beta(F(0)),
      result(result_),
      resultIndices(resultIndices_)
    {
    }

    virtual ~NonVoidOperation() {
    }

    virtual std::shared_ptr<Tensor<F>> getResult() {
      return result;
    }

    virtual const std::string &getResultIndices() {
      return resultIndices;
    }

  protected:
    F alpha, beta;
    std::shared_ptr<Tensor<F>> result;
    std::string resultIndices;

    friend class Tcc<F>;
  };
}

#endif

