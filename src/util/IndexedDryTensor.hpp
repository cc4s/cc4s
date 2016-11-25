/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef INDEXED_DRY_TENSOR_DEFINED
#define INDEXED_DRY_TENSOR_DEFINED

#include <util/DryTensorExpression.hpp>
#include <util/Log.hpp>
#include <string>

namespace cc4s {
  template <typename F=double>
  class DryTensor;

  template <typename F=double>
  class IndexedDryTensor: public cc4s::DryTensorExpression<F> {
  public:
    IndexedDryTensor(
      cc4s::DryTensor<F> const &tensor_, std::string const &indices_
    ): tensor(&tensor_), indices(indices_) {
    }
    virtual ~IndexedDryTensor() {
    }

    virtual void log() const {
      LOG(0, "TCC") << tensor->get_name() << "[" << indices << "]" << std::endl;
    }

    cc4s::DryTensor<F> &getTensor() const {
      return *tensor;
    }

    std::string const &getIndices() const {
      return indices;
    }

    /**
     * \brief Create another indexed tensor from this one having
     * different index names.
     **/
    IndexedDryTensor<F> const &operator [](std::string const &newIndices) {
      return IndexedDryTensor(tensor, newIndices);
    }

    /**
     * \brief Assigns the given right hand side expression to this
     * indexed tensor, returning this indexed tensor as result expression
     * for possible further operations.
     **/
    IndexedDryTensor<F> &operator =(
      cc4s::DryTensorExpression<F> const &T
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

