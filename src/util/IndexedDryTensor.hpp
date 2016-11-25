/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef INDEXED_DRY_TENSOR_DEFINED
#define INDEXED_DRY_TENSOR_DEFINED

#include <util/DryTensorTerm.hpp>
#include <util/Log.hpp>
#include <string>

namespace cc4s {
  template <typename F=double>
  class DryTensor;

  template <typename F=double>
  class IndexedDryTensor: public cc4s::DryTensorTerm<F> {
  public:
    IndexedDryTensor(
      cc4s::DryTensor<F> const &tensor_, std::string const &indices_
    ): tensor(&tensor_), indices(indices_) {
    }
    virtual ~IndexedDryTensor() {
    }

    virtual void log() const {
      LOG(0, "TCC") << tensor->location << "[" << indices << "]" << std::endl;
    }

    cc4s::DryTensor<F> &getTensor() const {
      return *tensor;
    }

    std::string const &getIndices() const {
      return indices;
    }

    cc4s::DryTensorTerm<F> &operator =(
      cc4s::DryTensorTerm<F> const &T
    ) {
      T.log();
      this->log();
      LOG(0,"TCC") << "summed" << std::endl;
      return *this;
    }

  protected:
    cc4s::DryTensor<F> const *tensor;
    std::string indices;
  };
}


#endif

