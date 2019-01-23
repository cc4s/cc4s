/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_INDEXING_DEFINED
#define TCC_INDEXING_DEFINED

#include <tcc/IndexedTensorExpression.hpp>

#include <util/SharedPointer.hpp>

#include <string>

namespace tcc {
  template <typename F, typename TE> ClosedTensorExpression;

  template <typename F, typename TE>
  class Indexing: public IndexedTensorExpression<F,TE> {
  public:
    /**
     * \brief Creates an expression with named indices from the closed tensor
     * source expression for further operations such as moves or contractions.
     * Not for direct invocation. Use the operator [] on
     * tcc::ClosedTensorExpression<F,TE> objects instead.
     **/
    Indexing(
      const PTR(ESC(ClosedTensorExpression<F,TE>)) &source_,
      const std::string &indices_,
      const typename Expression<TE>::ProtectedToken &
    ): source(source_), indices(indices_) {
    }

    virtual ~Indexing() {
    }

    /**
     * \brief Creates an expression with named indices from a closed tensor
     * source expression for further operations such as moves or contractions.
     * \param[in] source The closed tensor expression to operate on.
     * \param[in] indices The index character string where each character
     * specifies the index name of the respective dimension index in the
     * source expression.
     **/
    static PTR(ESC(Indexing<F,TE>)) create(
      const PTR(ESC(ClosedTensorExpression<F,TE>)) &source,
      const std::string &indices
    ) {
      return NEW(ESC(Indexing<F,TE>),
        source, indices, typename Expression<TE>::ProtectedToken()
      );
    }

    virtual PTR(Operation<TE>) compile(IndexCounts &indexCounts) {
//      return IndexingOperation<F,TE>::create(result, resultIndices);
    }

    virtual void countIndices(IndexCounts &indexCounts) {
      indexCounts.add(indices);
    }
  protected:
    PTR(ESC(ClosedTensorExpression<F,TE>)) source;
    std::string indices;
  };
}


#endif

