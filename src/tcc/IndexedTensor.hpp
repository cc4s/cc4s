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
  class IndexedTensor: public Expression<F> {
  public:
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

    std::shared_ptr<Tensor<F>> tensor;
    std::string indices;

  protected:
    /**
     * \brief Creates an expression with named indices from a stored
     * tensor for further operations such as contractions or assignments.
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

    friend class Tensor<F>;
  };

  /**
   * \brief Moves the given right hand side expression into the given
   * indexed tensor, returning the indexed tensor as result expression
   * for possible further operations.
   * Note that the operator = cannot be used since all expressions are
   * represented by shared pointers.
   **/
  template <typename Rhs>
  std::shared_ptr<Assignment<typename Rhs::FieldType>> operator <<=(
    const std::shared_ptr<IndexedTensor<typename Rhs::FieldType>> &lhs,
    const std::shared_ptr<Rhs> &rhs
  ) {
    return std::make_shared<Assignment<typename Rhs::FieldType>>(
      lhs, rhs
    );
  }

  template <typename Lhs, typename Rhs>
  std::shared_ptr<Assignment<typename Lhs::FieldType>> operator <<=(
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
    return std::shared_ptr<Assignment<typename Lhs::FieldType>>();
  }

  template <typename F>
  inline std::ostream &operator <<(
    std::ostream &stream, const IndexedTensor<F> &t
  ) {
    return stream << t.tensor->getName() << "[" << t.indices << "]";
  }
}


#endif

