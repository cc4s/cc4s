/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_TENSOR_RESULT_OPERATION_DEFINED
#define TCC_TENSOR_RESULT_OPERATION_DEFINED

#include <tcc/Operation.hpp>
#include <tcc/Costs.hpp>

#include <util/SharedPointer.hpp>
#include <util/Log.hpp>

#include <string>

namespace tcc {
  template <typename F> class Tensor;
  template <typename F> class Move;
  template <typename F> class Contraction;

  template <typename F>
  class TensorResultOperation: public Operation<F> {
  public:
    /**
     * \brief Creates a contraction operation contracting the results of
     * the left and the right sub-operations where the result is to be
     * stored in the specified result tensor.
     * Not intended for direct invocation. Use Tcc::compile(expression) to
     * generate operations.
     **/
    TensorResultOperation(
      const PTR(Tensor<F>) &result_,
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

    virtual ~TensorResultOperation() {
    }

    virtual PTR(Tensor<F>) getResult() {
      return result;
    }

    virtual const std::string &getResultIndices() {
      return resultIndices;
    }

  protected:
    F alpha, beta;
    PTR(Tensor<F>) result;
    std::string resultIndices;

    friend class Move<F>;
    friend class Contraction<F>;
  };
}

#endif

