#ifndef TCC_TENSOR_EXPRESSION_DEFINED
#define TCC_TENSOR_EXPRESSION_DEFINED

#include <tcc/Expression.hpp>

#include <tcc/TensorOperation.hpp>
#include <util/Exception.hpp>

namespace cc4s {
  template <typename F, typename TE>
  class TensorExpression: public Expression<TE> {
  public:
    typedef F FieldType;

    virtual PTR(ESC(TensorOperation<F,TE>)) lhsCompile(
      const PTR(ESC(TensorOperation<F,TE>)) &rhsOperation
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

