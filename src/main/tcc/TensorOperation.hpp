/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_TENSOR_OPERATION_DEFINED
#define TCC_TENSOR_OPERATION_DEFINED

#include <tcc/Operation.hpp>

#include <tcc/Costs.hpp>
#include <util/SharedPointer.hpp>
#include <util/Log.hpp>

namespace cc4s {
  template <typename F, typename TE> class Tensor;
  template <typename F, typename TE> class Slice;

  template <typename F, typename TE>
  class TensorOperation: public Operation<TE> {
  public:
    typedef F FieldType;

    TensorOperation(
      const PTR(ESC(Tensor<F,TE>)) &result_,
      const Costs &costs_,
      const std::string &file_, const size_t line_,
      const typename Operation<TE>::ProtectedToken &
    ):
      Operation<TE>(costs_, file_, line_),
      result(result_),
      alpha(F(1)), beta(F(0))
    {
    }

    void execute() override {
      // a tensor operation occurring as an atomic operation is a fetch
      // operation of the operand tensor, which result points to.
      // there is nothing more to do.
    }

    virtual PTR(ESC(Tensor<F,TE>)) getResult() const {
      return result;
    }

    template <typename G>
    bool isOlderThan(const Ptr<TensorOperation<G,TE>> &source) const {
      return
        // version is only relevant for replacing operations with beta==0
        beta != F(0) ||
        // equality means rhs & lhs tensor are the same, requiring update
        source->getResult()->getVersion() >= result->getVersion();
    }

    void updated() {
      result->updated();
    }

    std::string getName() const {
      return result->getName();
    }

  protected:
    PTR(ESC(Tensor<F,TE>)) result;
    F alpha, beta;

    friend class Tensor<F,TE>;
    friend class Slice<F,TE>;
  };
}

#endif

