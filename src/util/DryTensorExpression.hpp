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
}

// include all known expression types
#include <util/IndexedDryTensor.hpp>
#include <util/DryTensorContraction.hpp>

#endif

