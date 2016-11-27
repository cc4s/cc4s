/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_TENSOR_EXPRESSION_DEFINED
#define DRY_TENSOR_EXPRESSION_DEFINED

namespace cc4s {
  template <typename F=double>
  class DryTensorExpression {
  public:
    virtual ~DryTensorExpression() {
    }
    virtual void log() const = 0;
  protected:
  };

// forward class declaration of interdependent expression types
  template <typename F=double>
  class IndexedDryTensor;

  template <typename F=double>
  class DryTensorContraction;

  template <typename F=double>
  class DryTensorAssignment;
}


// include all known expression types
#include <tcc/IndexedDryTensor.hpp>
#include <tcc/DryTensorContraction.hpp>
#include <tcc/DryTensorAssignment.hpp>

#endif

