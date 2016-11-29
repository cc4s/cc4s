/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef INDEXED_TENSOR_DEFINED
#define INDEXED_TENSOR_DEFINED

#include <tcc/TensorExpression.hpp>
#include <util/StaticAssert.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <string>
#include <ostream>

namespace cc4s {
  template <typename F=double>
  class DryTensor;

  template <typename F>
  class IndexedTensor: public TensorExpression<F> {
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
    IndexedTensor(
      DryTensor<F> *tensor_, std::string const &indices_
    ): tensor(tensor_), indices(indices_) {
    }
    virtual ~IndexedTensor() {
    }

    virtual void log() const {
      LOG(0, "TCC") << tensor->get_name() << "[" << indices << "]" << std::endl;
    }

    virtual TensorOperation<F> *compile(std::string const &lhsIndices) {
      // not to be used
      throw new Exception("Operation required.");
    }

    /**
     * \brief Assigns the given right hand side expression to this
     * indexed tensor, returning this indexed tensor as result expression
     * for possible further operations.
     **/
    template <typename Rhs>
    TensorAssignment<typename Rhs::FieldType> &operator =(Rhs &rhs) {
      static_assert(
        TypeRelations<F, typename Rhs::FieldType>::Equals,
        "Assignment requires tensors of same type"
      );
      return *new TensorAssignment<typename Rhs::FieldType>(this, &rhs);
    }

    DryTensor<F> *tensor;
    std::string indices;
  };

  template <typename F>
  inline std::ostream &operator <<(
    std::ostream &stream, IndexedTensor<F> const &t
  ) {
    return stream << t.tensor->get_name() << "[" << t.indices << "]";
  }
}


#endif

