#ifndef TCC_INDEXED_TENSOR_OPERATION_DEFINED
#define TCC_INDEXED_TENSOR_OPERATION_DEFINED

#include <tcc/TensorOperation.hpp>

#include <tcc/Costs.hpp>
#include <util/SharedPointer.hpp>
#include <string>

namespace cc4s {
  template <typename F, typename TE> class Indexing;
  template <typename F, typename TE> class Move;
  template <typename F, typename TE> class Contraction;

  template <typename F, typename TE>
  class IndexedTensorOperation: public TensorOperation<F,TE> {
  public:
    typedef F FieldType;

    IndexedTensorOperation(
      const PTR(ESC(Tensor<F,TE>)) &result_,
      const char *resultIndices_,
      const Costs &costs_,
      const std::string &file_, const size_t line_,
      const typename Operation<TE>::ProtectedToken &
    ):
      TensorOperation<F,TE>(
        result_, costs_, file_, line_, typename Operation<TE>::ProtectedToken()
      ),
      resultIndices(resultIndices_)
    {
    }

    const std::string &getResultIndices() {
      return resultIndices;
    }

    std::string getName() const {
      return TensorOperation<F,TE>::getName() + "[" + resultIndices + "]";
    }

  protected:
    std::string resultIndices;

    friend class Indexing<F,TE>;
    friend class Move<F,TE>;
    friend class Contraction<F,TE>;
  };
}

#endif

