/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_TENSOR_RESULT_OPERATION_DEFINED
#define TCC_TENSOR_RESULT_OPERATION_DEFINED

#include <tcc/Operation.hpp>
#include <tcc/Costs.hpp>

#include <util/SharedPointer.hpp>

#include <string>

namespace tcc {
  template <typename F, typename TE> class Tensor;
  template <typename F, typename TE> class Move;
  template <typename F, typename TE> class Contraction;

  template <typename F, typename TE>
  class TensorResultOperation: public Operation<TE> {
  public:
    typedef F FieldType;

    /**
     * \brief Creates a contraction operation contracting the results of
     * the left and the right sub-operations where the result is to be
     * stored in the specified result tensor.
     * Not intended for direct invocation. Use Tcc::compile(expression) to
     * generate operations.
     **/
    TensorResultOperation(
      const PTR(ESC(Tensor<F,TE>)) &result_,
      const char *resultIndices_,
      const Costs &costs_,
      const typename Operation<TE>::ProtectedToken &
    ):
      Operation<TE>(costs_),
      result(result_),
      resultIndices(resultIndices_),
      alpha(F(1)), beta(F(0))
    {
    }

    virtual ~TensorResultOperation() {
    }

    virtual PTR(ESC(Tensor<F,TE>)) getResult() {
      return result;
    }

    virtual const std::string &getResultIndices() {
      return resultIndices;
    }

  protected:
    PTR(ESC(Tensor<F,TE>)) result;
    std::string resultIndices;
    F alpha, beta;

    friend class Move<F,TE>;
    friend class Contraction<F,TE>;
  };
}

#endif

