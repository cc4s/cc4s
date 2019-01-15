/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_INDEXED_TENSOR_DEFINED
#define TCC_INDEXED_TENSOR_DEFINED

#include <tcc/Expression.hpp>
#include <tcc/FetchOperation.hpp>
#include <util/SharedPointer.hpp>
#include <util/StaticAssert.hpp>

#include <string>
#include <ostream>

namespace tcc {
  template <typename F> class Tensor;
  template <typename F> class Move;

  template <typename F>
  class IndexedTensor: public Expression<F> {
  public:
    /**
     * \brief Creats an expression with named indices from a stored
     * tensor for further operations such as moves or contractions.
     * Not for direct invocation. Use tcc::IndexedTensor<F>::create
     * or the operator [] on tcc::Tensor<F> objects instead.
     **/
    IndexedTensor(
      const PTR(Tensor<F>) &tensor_, const std::string &indices_,
      const typename Expression<F>::ProtectedToken &protectedToken
    ): tensor(tensor_), indices(indices_) {
    }
    virtual ~IndexedTensor() {
    }

    /**
     * \brief Creates an expression with named indices from a stored
     * tensor for further operations such as moves or contractions.
     * \param[in] tensor The stored tensor to be operated on.
     * \param[in] indices The index character string where each character
     * specifies the index name of the respective dimension index in this
     * expression.
     **/
    static PTR(IndexedTensor<F>) create(
      const PTR(Tensor<F>) &tensor, const std::string &indices
    ) {
      return NEW(IndexedTensor<F>,
        tensor, indices, typename Expression<F>::ProtectedToken()
      );
    }

    virtual PTR(Operation<F>) compile(IndexCounts &indexCounts) {
      return FetchOperation<F>::create(DYNAMIC_PTR_CAST(IndexedTensor<F>,THIS));
    }

    virtual void countIndices(IndexCounts &indexCounts) {
      indexCounts.add(indices);
    }

    PTR(Tensor<F>) getTensor() {
      return tensor;
    }

    std::string getIndices() {
      return indices;
    }

  protected:
    PTR(Tensor<F>) tensor;
    std::string indices;

    friend class Move<F>;
  };
}


#endif

