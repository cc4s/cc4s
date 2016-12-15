/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_INDEXED_TENSOR_DEFINED
#define TCC_INDEXED_TENSOR_DEFINED

#include <tcc/Expression.hpp>
#include <tcc/Move.hpp>
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
  class IndexedTensor: public Expression<F> {
  public:
    /**
     * \brief Creats an expression with named indices from a stored
     * tensor for further operations such as moves or contractions.
     * Not for direct invocation. Use tcc::IndexedTensor<F>::create
     * or the operator [] on tcc::Tensor<F> objects instead.
     **/
    IndexedTensor(
      const std::shared_ptr<Tensor<F>> &tensor_, const std::string &indices_,
      const typename Expression<F>::ProtectedToken &protectedToken
    ): tensor(tensor_), indices(indices_) {
    }
    virtual ~IndexedTensor() {
    }

    virtual std::shared_ptr<Operation<F>> compile(
      std::string const &lhsIndices
    ) {
      // not to be used
      throw new EXCEPTION("Operation on IndexedTensor required.");
    }

    /**
     * \brief Creates an expression with named indices from a stored
     * tensor for further operations such as moves or contractions.
     * \param[in] tensor The stored tensor to be operated on.
     * \param[in] indices The index character string where each character
     * specifies the index name of the respective dimension index in this
     * expression.
     **/
    static std::shared_ptr<IndexedTensor<F>> create(
      const std::shared_ptr<Tensor<F>> &tensor, const std::string &indices
    ) {
      return std::make_shared<IndexedTensor<F>>(
        tensor, indices, typename Expression<F>::ProtectedToken()
      );
    }

    std::shared_ptr<Tensor<F>> tensor;
    std::string indices;

  protected:
    template <typename Rhs>
    friend std::shared_ptr<Move<typename Rhs::FieldType>> operator <<=(
      const std::shared_ptr<IndexedTensor<typename Rhs::FieldType>> &lhs,
      const std::shared_ptr<Rhs> &rhs
    );
  };

  /**
   * \brief Moves the given right hand side expression into the given
   * indexed tensor, returning the indexed tensor as result expression
   * for possible further operations.
   * Note that the operator = cannot be used since all expressions are
   * represented by shared pointers.
   **/
  template <typename Rhs>
  inline std::shared_ptr<Move<typename Rhs::FieldType>> operator <<=(
    const std::shared_ptr<IndexedTensor<typename Rhs::FieldType>> &lhs,
    const std::shared_ptr<Rhs> &rhs
  ) {
    return std::make_shared<Move<typename Rhs::FieldType>>(
      lhs, rhs,
      typename Expression<typename Rhs::FieldType>::ProtectedToken()
    );
  }

  template <typename Lhs, typename Rhs>
  inline std::shared_ptr<Move<typename Lhs::FieldType>> operator <<=(
    const std::shared_ptr<Lhs> &, const std::shared_ptr<Rhs> &rhs
  ) {
    static_assert(
      cc4s::TypeRelations<
        typename Lhs::FieldType, typename Rhs::FieldType
      >::Equals,
      "Move operations requires tensors of same type."
    );
    static_assert(
      cc4s::StaticAssert<Lhs>::False,
      "Only indexed tensors may be used as the left hand side of a move operation."
    );
    return std::shared_ptr<Move<typename Lhs::FieldType>>();
  }

  template <typename F>
  inline std::ostream &operator <<(
    std::ostream &stream, const IndexedTensor<F> &t
  ) {
    return stream << t.tensor->getName() << "[" << t.indices << "]";
  }
}


#endif

