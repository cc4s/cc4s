/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_TENSOR_EXPRESSION_DEFINED
#define TCC_TENSOR_EXPRESSION_DEFINED

#include <tcc/Expression.hpp>
#include <tcc/TensorOperation.hpp>

namespace tcc {
  template <typename F, typename TE>
  class TensorExpression: public Expression<TE> {
  public:
    typedef F FieldType;

    virtual PTR(ESC(TensorOperation<F,TE>)) lhsCompile(
      const PTR(ESC(TensorOperation<F,TE>)) &rhsOperation
    ) {
      throw new EXCEPTION(
        "Writable expression expected on the "
        "left-hand-side of " + rhsOperation->getResult()->getName()
      );
    }
  };
}

#endif

