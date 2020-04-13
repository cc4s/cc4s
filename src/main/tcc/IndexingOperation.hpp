/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_INDEXING_OPERATION_DEFINED
#define TCC_INDEXING_OPERATION_DEFINED

#include <tcc/IndexedTensorOperation.hpp>

#include <tcc/Costs.hpp>
#include <util/SharedPointer.hpp>
#include <string>

namespace cc4s {
  template <typename F, typename TE> class Indexing;

  template <typename F, typename TE>
  class IndexingOperation: public IndexedTensorOperation<F,TE> {
  public:
    typedef F FieldType;

    IndexingOperation(
      const PTR(ESC(TensorOperation<F,TE>)) &source_,
      const PTR(ESC(Tensor<F,TE>)) &result_,
      const char *resultIndices_,
      const Costs &costs_,
      const typename Operation<TE>::ProtectedToken &
    ):
      IndexedTensorOperation<F,TE>(
        result_, resultIndices_, costs_,
        typename Operation<TE>::ProtectedToken()
      ),
      source(source_)
    {
    }

    virtual ~IndexingOperation() {
    }

    virtual void execute() {
      source->execute();
    }

    virtual operator std::string () const {
      return std::string(*source) + "[" + this->resultIndices + "]";
    }

  protected:
    PTR(ESC(TensorOperation<F,TE>)) source;

    static PTR(ESC(IndexingOperation<F,TE>)) create(
      const PTR(ESC(TensorOperation<F,TE>)) &source,
      const PTR(ESC(Tensor<F,TE>)) &result,
      const char *resultIndices,
      const Costs &costs
    ) {
      return NEW(ESC(IndexingOperation<F,TE>),
        source, result, resultIndices, costs,
        typename Operation<TE>::ProtectedToken()
      );
    }

    friend class Indexing<F,TE>;
  };
}

#endif

