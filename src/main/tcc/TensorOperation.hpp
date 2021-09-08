#ifndef TCC_TENSOR_OPERATION_DEFINED
#define TCC_TENSOR_OPERATION_DEFINED

#include <tcc/Operation.hpp>
#include <tcc/Costs.hpp>
#include <math/Real.hpp>
#include <math/Complex.hpp>
#include <util/SharedPointer.hpp>
#include <util/Log.hpp>

namespace cc4s {
  template <typename F, typename TE> class Tensor;
  template <typename F, typename TE> class Slice;
  template <typename F> class TensorOperationTraits;

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
      // FIXME: discard versioning based on isOlderThan
      return true;
    }

    size_t getLatestSourceVersion() override {
      return result->getVersion();
    }

    void updated() {
      result->updated();
    }

    std::string getName() const {
      return result->getName();
    }

    PTR(ESC(Tensor<F,TE>)) getResult() {
      return result;
    }


  protected:
    PTR(ESC(Tensor<F,TE>)) result;
    F alpha, beta;

    void accountFlops() {
      Operation<TE>::flops += this->costs.additionsCount *
        TensorOperationTraits<F>::getFlopsPerAddition();
      Operation<TE>::flops += this->costs.multiplicationsCount *
        TensorOperationTraits<F>::getFlopsPerMultiplication();
    }

    friend class Tensor<F,TE>;
    friend class Slice<F,TE>;
  };

  template <>
  class TensorOperationTraits<Real<64>> {
  public:
    static size_t getFlopsPerAddition() { return 1; }
    static size_t getFlopsPerMultiplication() { return 1; }
  };
  template <>
  class TensorOperationTraits<Complex<64>> {
  public:
    static size_t getFlopsPerAddition() { return 2; }
    static size_t getFlopsPerMultiplication() { return 6; }
  };
}

#endif

