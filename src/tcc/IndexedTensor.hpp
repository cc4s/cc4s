/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_INDEXED_TENSOR_DEFINED
#define TCC_INDEXED_TENSOR_DEFINED

#include <tcc/Expression.hpp>
#include <util/StaticAssert.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>

#include <string>
#include <ostream>
#include <memory>

namespace tcc {
  template <typename F=double>
  class Tensor;

  template <typename F>
  class IndexedTensor:
    public Expression<F>,
    public std::enable_shared_from_this<IndexedTensor<F>>
  {
  public:
    /**
     * \brief Creates an expression with named indices from a stored
     * tensor for further operations such as contractions or assignments.
     * \param[in] tensor_ The stored tensor to be operated on.
     * \param[in] indices_ The index character string where each character
     * specifies the index name of the respective dimension index in this
     * expression.
     **/
    IndexedTensor(
      const std::shared_ptr<Tensor<F>> &tensor_, const std::string &indices_
    ): tensor(tensor_), indices(indices_) {
    }
    virtual ~IndexedTensor() {
    }

    virtual std::shared_ptr<Operation<F>> compile(
      std::string const &lhsIndices
    ) {
      // not to be used
      throw new EXCEPTION("Operation required.");
    }

    /**
     * \brief Assigns the given right hand side expression to this
     * indexed tensor, returning this indexed tensor as result expression
     * for possible further operations.
     **/
    template <typename Rhs>
    Assignment<typename Rhs::FieldType> &operator =(
      const std::shared_ptr<Rhs> &rhs
    ) {
      static_assert(
        cc4s::TypeRelations<F, typename Rhs::FieldType>::Equals,
        "Assignment requires tensors of same type"
      );
      return std::make_shared<Assignment<typename Rhs::FieldType>>(
        this->shared_from_this(), rhs
      );
    }

    std::shared_ptr<Tensor<F>> tensor;
    std::string indices;
  };

  template <typename F>
  inline std::ostream &operator <<(
    std::ostream &stream, const IndexedTensor<F> &t
  ) {
    return stream << t.tensor->get_name() << "[" << t.indices << "]";
  }
}


#endif

