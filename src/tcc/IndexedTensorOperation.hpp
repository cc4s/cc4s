/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_INDEXED_TENSOR_OPERATION_DEFINED
#define TCC_INDEXED_TENSOR_OPERATION_DEFINED

#include <tcc/TensorOperation.hpp>
#include <tcc/Costs.hpp>

#include <util/SharedPointer.hpp>

#include <string>

namespace tcc {
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
      const typename Operation<TE>::ProtectedToken &
    ):
      TensorOperation<F,TE>(
        result_, costs_, typename Operation<TE>::ProtectedToken()
      ),
      resultIndices(resultIndices_)
    {
    }

    virtual ~IndexedTensorOperation() {
    }

    virtual const std::string &getResultIndices() {
      return resultIndices;
    }

  protected:
    std::string resultIndices;

    static PTR(ESC(IndexedTensorOperation<F,TE>)) create(
      const PTR(ESC(Tensor<F,TE>)) &result,
      const char *resultIndices,
      const Costs &costs
    ) {
      return NEW(ESC(IndexedTensorOperation<F,TE>),
        result, resultIndices, costs, typename Operation<TE>::ProtectedToken()
      );
    }

    friend class Indexing<F,TE>;
    friend class Move<F,TE>;
    friend class Contraction<F,TE>;
  };
}

#endif

