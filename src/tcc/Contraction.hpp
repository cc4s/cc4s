/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TCC_CONTRACTION_DEFINED
#define TCC_CONTRACTION_DEFINED

#include <tcc/Expression.hpp>
#include <util/StaticAssert.hpp>

#include <vector>
#include <memory>

namespace tcc {
  template <typename F>
  class IndexedTensor;

  template <typename F>
  class Contraction: public Expression<F> {
  public:
    /**
     * \brief Creates a contraction expression of the two given tensor
     * expressions A and B.
     **/
    template <typename Lhs, typename Rhs>
    static std::shared_ptr<Contraction<typename Lhs::FieldType>> create(
      const std::shared_ptr<Lhs> &A, const std::shared_ptr<Rhs> &B
    ) {
      static_assert(
        cc4s::TypeRelations<
          typename Lhs::FieldType, typename Rhs::FieldType
        >::EQUALS,
        "Only tensors of the same type can be contracted."
      );
      auto contraction(
        std::make_shared<Contraction<typename Lhs::FieldType>>(
          A, B,
          typename Expression<typename Lhs::FieldType>::ProtectedToken()
        )
      );
      A->parent = contraction;
      B->parent = contraction;
      return contraction;
    }

    /**
     * \brief Creates a contraction expression of a tensor expression A
     * and a scalar alpha
     **/
    template <typename S, typename Lhs>
    static std::shared_ptr<Contraction<typename Lhs::FieldType>> create(
      const S alpha, const std::shared_ptr<Lhs> &A
    ) {
      static_assert(
        cc4s::TypeRelations<S, typename Lhs::FieldType>::CASTABLE_TO,
        "The type of the scalar must be convertible to the tensor type."
      );
      auto contraction(
        std::make_shared<Contraction<typename Lhs::FieldType>>(
          alpha, A,
          typename Expression<typename Lhs::FieldType>::ProtectedToken()
        )
      );
      A->parent = contraction;
      return contraction;
    }


    /**
     * \brief Flattening constructor given two contractions.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const std::shared_ptr<Contraction<F>> &lhs,
      const std::shared_ptr<Contraction<F>> &rhs,
      const typename Expression<F>::ProtectedToken &
    ): alpha(lhs->alpha * rhs->alpha), factors(lhs->factors) {
      factors.insert(factors.end(), rhs->factors.begin(), rhs->factors.end());
    }
    /**
     * \brief Flattening constructor given a contraction on the left hand
     * side and another expression on the right hand side.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const std::shared_ptr<Contraction<F>> &lhs,
      const std::shared_ptr<IndexedTensor<F>> &rhs,
      const typename Expression<F>::ProtectedToken &
    ): alpha(lhs->alpha), factors(lhs->factors) {
      factors.push_back(rhs);
    }
    /**
     * \brief Flattening constructor given a contraction on the right hand
     * side and another expression on the left hand side.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const std::shared_ptr<IndexedTensor<F>> &lhs,
      const std::shared_ptr<Contraction<F>> &rhs,
      const typename Expression<F>::ProtectedToken &
    ): alpha(rhs->alpha), factors(rhs->factors) {
      factors.push_back(lhs);
    }
    /**
     * \brief Constructor given two indexed tensors.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const std::shared_ptr<IndexedTensor<F>> &lhs,
      const std::shared_ptr<IndexedTensor<F>> &rhs,
      const typename Expression<F>::ProtectedToken &
    ): alpha(F(1)) {
      factors.push_back(lhs);
      factors.push_back(rhs);
    }
    /**
     * \brief Constructor given two general expressions.
     * This is currently not supported.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const std::shared_ptr<Expression<F>> &lhs,
      const std::shared_ptr<Expression<F>> &rhs,
      const typename Expression<F>::ProtectedToken &
    ) {
      static_assert(
        cc4s::StaticAssert<F>::FALSE,
        "Only contractions of contractions or tensors supported."
      );
    }
    /**
     * \brief Constructor given an indexed tensor and a scalar.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const F alpha_,
      const std::shared_ptr<IndexedTensor<F>> &lhs,
      const typename Expression<F>::ProtectedToken &
    ): alpha(alpha_) {
      factors.push_back(lhs);
    }
    /**
     * \brief Falttening constructor given a contraction and a scalar.
     * Not intended for direct invocation, create contractions
     * using the static create method.
     **/
    Contraction(
      const F alpha_,
      const std::shared_ptr<Contraction<F>> &lhs,
      const typename Expression<F>::ProtectedToken &
    ): alpha(lhs->alpha * alpha_), factors(lhs->factors) {
    }

    virtual ~Contraction() {
    }

    F alpha;
    std::vector<std::shared_ptr<IndexedTensor<F>>> factors;
  };

  /**
   * \brief Creates a contraction expression of the two given tensor
   * expressions A and B using the multiplication operator *.
   **/
  template <typename Lhs, typename Rhs>
  inline std::shared_ptr<Contraction<typename Lhs::FieldType>> operator *(
    const std::shared_ptr<Lhs> &A, const std::shared_ptr<Rhs> &B
  ) {
    return Contraction<typename Lhs::FieldType>::create(A, B);
  }

  /**
   * \brief Creates a contraction expression of a given tensor
   * expressions A and a scalar alpha using the multiplication operator *.
   **/
  template <typename Lhs, typename S>
  inline std::shared_ptr<Contraction<typename Lhs::FieldType>> operator *(
    const std::shared_ptr<Lhs> &A, const S alpha
  ) {
    return Contraction<typename Lhs::FieldType>::create(alpha, A);
  }

  /**
   * \brief Creates a contraction expression of a given tensor
   * expressions A and a scalar alpha using the multiplication operator *.
   **/
  template <typename S, typename Rhs>
  inline std::shared_ptr<Contraction<typename Rhs::FieldType>> operator *(
    const S alpha, const std::shared_ptr<Rhs> &A
  ) {
    return Contraction<typename Rhs::FieldType>::create(alpha, A);
  }
}

#endif

