/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_TENSOR_TERM_DEFINED
#define DRY_TENSOR_TERM_DEFINED

namespace cc4s {
  template <typename F=double>
  class DryTensorTerm {
  public:
    virtual ~DryTensorTerm() {
    }
    virtual void log() const = 0;
  protected:
  };
}

#include <util/DryTensorContraction.hpp>
#include <util/IndexedDryTensor.hpp>

#endif

