/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_CLOSED_TENSOR_EXPRESSION_DEFINED
#define TCC_CLOSED_TENSOR_EXPRESSION_DEFINED

#include <tcc/TensorExpression.hpp>

#include <tcc/Indexing.hpp>
#include <tcc/Slice.hpp>
#include <util/SharedPointer.hpp>
#include <initializer_list>

namespace cc4s {
  /**
   * \brief Closed tensor expressions are tensor valued expressions
   * whose dimensions are referred to by their slot number rather than
   * by a specific index label. Tensors and slices of tensors are
   * examples.
   **/
  template <typename F, typename TE>
  class ClosedTensorExpression: public TensorExpression<F,TE> {
  public:
    virtual void countIndices(Scope &scope) {
      // a closed tensor expression has no indices
    }

    /**
     * \brief Specify named indices of this tensor to be used in a
     * tensor expression. Indexed tensors are atomic types of tensor
     * expressions.
     **/
    PTR(ESC(Indexing<F,TE>)) operator[](
      const std::string &indices
    ) {
      return Indexing<F,TE>::create(
        THIS(ESC(ClosedTensorExpression<F,TE>)), indices
      );
    }

    PTR(ESC(Slice<F,TE>)) operator()(
      std::initializer_list<size_t> min,
      std::initializer_list<size_t> max
    ) {
      return Slice<F,TE>::create(
        THIS(ESC(ClosedTensorExpression<F,TE>)),
        min, max
      );
    }
  };
}

#endif

