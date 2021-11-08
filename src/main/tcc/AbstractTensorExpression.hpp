/*Copyright (c) 2019, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_ABSTRACT_TENSOR_EXPRESSION_DEFINED
#define TCC_ABSTRACT_TENSOR_EXPRESSION_DEFINED

#include <tcc/Expression.hpp>

#include <tcc/TensorOperation.hpp>
#include <util/Exception.hpp>

namespace cc4s {
  template <typename F, typename TE>
  class AbstractTensorExpression: public Expression<TE> {
  public:
    typedef F FieldType;

    virtual Ptr<TensorOperation<F,TE>> lhsCompile(
      const Ptr<TensorOperation<F,TE>> &rhsOperation
    ) {
      throw New<Exception>(
        "Writable expression expected on the "
        "left-hand-side of " + rhsOperation->getResult()->getName(),
        SourceLocation(rhsOperation->file, rhsOperation->line)
      );
    }
  };
}

#endif

