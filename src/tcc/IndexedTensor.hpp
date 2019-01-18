/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_INDEXED_TENSOR_DEFINED
#define TCC_INDEXED_TENSOR_DEFINED

#include <tcc/TensorResultExpression.hpp>
#include <tcc/FetchOperation.hpp>
#include <util/SharedPointer.hpp>
#include <util/StaticAssert.hpp>

#include <string>
#include <ostream>

namespace tcc {
  template <typename F, typename TE> class Tensor;
  template <typename F, typename TE> class Move;

  template <typename F, typename TE>
  class IndexedTensor: public TensorResultExpression<F,TE> {
  public:
    /**
     * \brief Creats an expression with named indices from a stored
     * tensor for further operations such as moves or contractions.
     * Not for direct invocation. Use tcc::IndexedTensor<F,TE>::create
     * or the operator [] on tcc::Tensor<F> objects instead.
     **/
    IndexedTensor(
      const PTR(ESC(Tensor<F,TE>)) &tensor_,
      const std::string &indices_,
      const typename Expression<TE>::ProtectedToken &protectedToken
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
    static PTR(ESC(IndexedTensor<F,TE>)) create(
      const PTR(ESC(Tensor<F,TE>)) &tensor,
      const std::string &indices
    ) {
      return NEW(ESC(IndexedTensor<F,TE>),
        tensor, indices, typename Expression<TE>::ProtectedToken()
      );
    }

    virtual PTR(Operation<TE>) compile(IndexCounts &indexCounts) {
      return FetchOperation<F,TE>::create(tensor, indices);
    }

    virtual void countIndices(IndexCounts &indexCounts) {
      indexCounts.add(indices);
    }

    PTR(ESC(Tensor<F,TE>)) getTensor() {
      return tensor;
    }

    std::string getIndices() {
      return indices;
    }

  protected:
    PTR(ESC(Tensor<F,TE>)) tensor;
    std::string indices;

    friend class Move<F,TE>;
  };
}


#endif

