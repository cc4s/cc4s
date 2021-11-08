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
      const Ptr<TensorOperation<F,TE>> &source_,
      const Ptr<Tensor<F,TE>> &result_,
      const char *resultIndices_,
      const Costs &costs_,
      const std::string &file_, const size_t line_,
      const typename Operation<TE>::ProtectedToken &
    ):
      IndexedTensorOperation<F,TE>(
        result_, resultIndices_, costs_,
        file_, line_, typename Operation<TE>::ProtectedToken()
      ),
      source(source_)
    {
    }

    void execute() override {
      source->execute();
    }

    operator std::string () const override {
      return std::string(*source) + "[" + this->resultIndices + "]";
    }

  protected:
    Ptr<TensorOperation<F,TE>> source;

    static Ptr<IndexingOperation<F,TE>> create(
      const Ptr<TensorOperation<F,TE>> &source,
      const Ptr<Tensor<F,TE>> &result,
      const char *resultIndices,
      const Costs &costs,
      const Scope &scope
    ) {
      return New<IndexingOperation<F,TE>>(
        source, result, resultIndices, costs,
        scope.file, scope.line, typename Operation<TE>::ProtectedToken()
      );
    }

    friend class Indexing<F,TE>;
  };
}

#endif

