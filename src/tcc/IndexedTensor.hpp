/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_INDEXED_TENSOR_DEFINED
#define TCC_INDEXED_TENSOR_DEFINED

#include <tcc/Expression.hpp>
#include <tcc/Move.hpp>
#include <util/SharedPointer.hpp>
#include <util/StaticAssert.hpp>

#include <string>
#include <ostream>

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
      const PTR(Tensor<F>) &tensor_, const std::string &indices_,
      const typename Expression<F>::ProtectedToken &protectedToken
    ): tensor(tensor_), indices(indices_) {
    }
    virtual ~IndexedTensor() {
    }

    /**
     * \brief Creates an expression with named indices from a stored
     * tensor for further operations such as moves or contractions.
     * \param[in] tensor The stored tensor to be operated on.
     * \param[in] indices The index character string where each character
     * specifies the index name of the respective dimension index in this
     * expression.
     **/
    static PTR(IndexedTensor<F>) create(
      const PTR(Tensor<F>) &tensor, const std::string &indices
    ) {
      return NEW(IndexedTensor<F>,
        tensor, indices, typename Expression<F>::ProtectedToken()
      );
    }

    PTR(Tensor<F>) tensor;
    std::string indices;
  };

  /**
   * \brief Moves the given right hand side expression into the given
   * indexed tensor, returning the indexed tensor as result expression
   * for possible further operations.
   * Note that the operator = cannot be used since all expressions are
   * represented by shared pointers.
   **/
  template <typename Rhs>
  inline PTR(Move<typename Rhs::FieldType>) operator <<=(
    const PTR(IndexedTensor<typename Rhs::FieldType>) &lhs,
    const PTR(Rhs) &rhs
  ) {
    return Move<typename Rhs::FieldType>::create(lhs, rhs, 0);
  }

  template <typename Lhs, typename Rhs>
  inline PTR(Move<typename Lhs::FieldType>) operator <<=(
    const PTR(Lhs) &, const PTR(Rhs) &rhs
  ) {
    static_assert(
      cc4s::TypeRelations<
        typename Lhs::FieldType, typename Rhs::FieldType
      >::EQUALS,
      "Move operations requires tensors of same type."
    );
    static_assert(
      cc4s::StaticAssert<Lhs>::FALSE,
      "Only indexed tensors may be used as the left hand side of a move operation."
    );
    return PTR(Move<typename Lhs::FieldType>)();
  }

  template <typename F>
  inline std::ostream &operator <<(
    std::ostream &stream, const IndexedTensor<F> &t
  ) {
    return stream << t.tensor->getName() << "[" << t.indices << "]";
  }
}


#endif

