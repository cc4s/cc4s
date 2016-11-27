/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef INDEXED_DRY_TENSOR_DEFINED
#define INDEXED_DRY_TENSOR_DEFINED

#include <tcc/DryTensorExpression.hpp>
#include <util/Log.hpp>
#include <string>

namespace cc4s {
  template <typename F=double>
  class DryTensor;

  template <typename F>
  class IndexedDryTensor: public DryTensorExpression<F> {
  public:
    /**
     * \brief Creates an expression with named indices from a stored
     * tensor further operations such as contractions or index permutations.
     * \param[in] tensor_ The stored tensor to be operated on.
     * The tensor is an idependent entity and will not be destroyed upon
     * destruction of this expression.
     * \param[in] indices_ The index character string where each character
     * specifies the index name of the respective dimension index in this
     * expression.
     **/
    IndexedDryTensor(
      DryTensor<F> *tensor_, std::string const &indices_
    ): tensor(tensor_), indices(indices_) {
    }
    virtual ~IndexedDryTensor() {
    }

    virtual void log() const {
      LOG(0, "TCC") << tensor->get_name() << "[" << indices << "]" << std::endl;
    }

    /**
     * \brief Assigns the given right hand side expression to this
     * indexed tensor, returning this indexed tensor as result expression
     * for possible further operations.
     **/
    template <typename Rhs>
    DryTensorAssignment<typename Rhs::FieldType> &operator =(Rhs &rhs) {
      return *new DryTensorAssignment<typename Rhs::FieldType>(this, &rhs);
    }

    DryTensor<F> *tensor;
    std::string indices;
  };
}


#endif

