/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_INDEXED_TENSOR_EXPRESSION_DEFINED
#define TCC_INDEXED_TENSOR_EXPRESSION_DEFINED

#include <tcc/TensorExpression.hpp>

namespace cc4s {
  /**
   * \brief Indexed tensor expressions are tensor valued expressions
   * whose dimensions are referred to by a specific index labels
   * for indentifying related indices in other indexed tensors.
   * Results of contractions are examples.
   * Indexing and closing converts between closed and indexed tensor
   * expressions. The transposition of a matrix tensor shared pointer A is
   * for instance done by: (*A)["ij"]^"ji"
   **/
  template <typename F, typename TE>
  class IndexedTensorExpression: public TensorExpression<F,TE> {
  public:
  };
}

#endif

