/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_INDEXING_DEFINED
#define TCC_INDEXING_DEFINED

#include <tcc/IndexedTensorExpression.hpp>

#include <tcc/IndexingOperation.hpp>
#include <util/SharedPointer.hpp>
#include <string>

namespace tcc {
  template <typename F, typename TE> class ClosedTensorExpression;

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

    virtual PTR(Operation<TE>) compile(Scope &scope) {
      auto sourceOperation(
        DYNAMIC_PTR_CAST(
          ESC(TensorOperation<F,TE>), source->compile()
        )
      );
      return IndexingOperation<F,TE>::create(
        sourceOperation,
        sourceOperation->getResult(), indices.c_str(), sourceOperation->costs
      );
    }

    // keep other overloads visible
    using Expression<TE>::compile;

    virtual PTR(ESC(TensorOperation<F,TE>)) lhsCompile(
      const PTR(ESC(TensorOperation<F,TE>)) &rhsOperation
    )  {
      auto indexedRhs(
        DYNAMIC_PTR_CAST(ESC(IndexedTensorOperation<F,TE>), rhsOperation)
      );
      if (indexedRhs) {
        // change index order of indexed operation result according to
        // index order provided by this indexing expression
        transposeResult(indexedRhs);
        // continue lhs compiling of contained closed tensor expression
        return source->lhsCompile(indexedRhs);
      }
      throw new EXCEPTION(
        "Expecting indexed expression ([\"...\"]) on the left-hand-side of "
        "an assignment if the right-hand-side carries indices."
      );
    }

    virtual void countIndices(Scope &scope) {
      scope.add(indices);
    }

    PTR(ESC(ClosedTensorExpression<F,TE>)) source;
    std::string indices;

  protected:
    // transposes the result of the indexed rhs tensor to match the index
    // order of this lhs indexed expression
    void transposeResult(
      const PTR(ESC(IndexedTensorOperation<F,TE>)) &rhs
    ) {
      if (indices.length() != rhs->getResultIndices().length()) {
        throw new EXCEPTION(
          "Number of indices of left-hand-side expression "
          "must match the number of indices of the right-hand-side result."
        );
      }
      // determine shape of lhs
      size_t lenOfIndex[std::numeric_limits<char>::max()+1];
      // enter which indices correspond to which length on the rhs
      for (unsigned int i(0); i < indices.length(); ++i) {
        lenOfIndex[static_cast<unsigned int>(rhs->getResultIndices()[i])]=
          rhs->getResult()->lens[i];
      }
      // transpose rhs result
      for (unsigned int i(0); i < indices.length(); ++i) {
        rhs->getResult()->lens[i] = lenOfIndex[
          static_cast<unsigned int>(indices[i])
        ];
      }
      rhs->resultIndices = indices;
    }
  };
}


#endif

