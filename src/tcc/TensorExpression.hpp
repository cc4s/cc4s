/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_TENSOR_EXPRESSION_DEFINED
#define TCC_TENSOR_EXPRESSION_DEFINED

#include <tcc/Expression.hpp>

namespace tcc {
  template <typename F, typename TE>
  class TensorExpression: public Expression<TE> {
  public:
    typedef F FieldType;

    virtual std::vector<size_t> getLens() = 0;
  };
}

#endif

