/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_TENSOR_EXPRESSION_DEFINED
#define DRY_TENSOR_EXPRESSION_DEFINED

#include <util/StaticAssert.hpp>

namespace cc4s {
  // forward class declaration of interdependent expression types
  template <typename F=double>
  class DryTensorExpression;

  template <typename F=double>
  class IndexedDryTensor;

  template <typename F=double>
  class DryTensorContraction;

  template <typename F=double>
  class DryTensorAssignment;

  template <typename F>
  class DryTensorExpression {
  public:
    typedef F FieldType;
    virtual ~DryTensorExpression() {
    }
    virtual void log() const = 0;

    template <typename Rhs>
    DryTensorAssignment<F> &operator =(Rhs &rhs) {
      static_assert(
        StaticAssert<F>::False,
        "Only indexed tensors may be used as the left hand side of an assignment."
      );
    }
  protected:
  };
}


// include all known expression types
#include <tcc/IndexedDryTensor.hpp>
#include <tcc/DryTensorContraction.hpp>
#include <tcc/DryTensorAssignment.hpp>

#endif

