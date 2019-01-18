/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_TENSOR_RESULT_EXPRESSION_DEFINED
#define TCC_TENSOR_RESULT_EXPRESSION_DEFINED

#include <tcc/Expression.hpp>

namespace tcc {
  template <typename F, typename TE>
  class TensorResultExpression: public Expression<TE> {
  public:
    typedef F FieldType;
  };
}

#endif

