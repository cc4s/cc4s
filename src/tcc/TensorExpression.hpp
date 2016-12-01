/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_EXPRESSION_DEFINED
#define TENSOR_EXPRESSION_DEFINED

#include <util/StaticAssert.hpp>

#include <memory>
using std::shared_ptr;

// TODO: use shared pointers for entire tcc

namespace cc4s {
  // forward class declaration of interdependent expression types
  template <typename F=double>
  class TensorExpression;

  template <typename F=double>
  class IndexedTensor;

  template <typename F=double>
  class TensorContraction;

  template <typename F=double>
  class TensorAssignment;

  template <typename F=double>
  class TensorOperation;

  template <typename F>
  class TensorExpression {
  public:
    typedef F FieldType;
    virtual ~TensorExpression() {
    }

    template <typename Rhs>
    TensorAssignment<F> &operator =(Rhs &rhs) {
      static_assert(
        StaticAssert<F>::False,
        "Only indexed tensors may be used as the left hand side of an assignment."
      );
    }

    /**
     * \brief Compile this tensor expression into a TensorOperation.
     **/
    virtual shared_ptr<TensorOperation<F>> compile(
      std::string const &lhsIndices
    ) = 0;
  };
}

#endif

