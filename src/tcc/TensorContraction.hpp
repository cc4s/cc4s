/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TENSOR_CONTRACTION_DEFINED
#define TENSOR_CONTRACTION_DEFINED

#include <tcc/TensorExpression.hpp>
#include <util/StaticAssert.hpp>
#include <util/Log.hpp>
#include <vector>

namespace cc4s {
  template <typename F>
  class TensorContraction: public TensorExpression<F> {
  public:
    /**
     * \brief Flattening constructor given two contractions.
     **/
    TensorContraction(
      TensorContraction<F> *lhs, TensorContraction<F> *rhs
    ): factors(lhs.factors) {
      factors.insert(factors.end(), rhs->factors.begin(), rhs->factors.end());
      // the factors from both lhs and rhs contractions are now contained here
      lhs->factors.clear(); rhs->factors.clear();
      delete lhs, rhs;
    }
    /**
     * \brief Flattening constructor given a contraction on the left hand
     * side and another expression on the right hand side.
     **/
    TensorContraction(
      TensorContraction<F> *lhs, IndexedTensor<F> *rhs
    ): factors(lhs->factors) {
      factors.push_back(rhs);
      // the factors from the lhs expression are now contained here
      lhs->factors.clear();
      delete lhs;
    }
    /**
     * \brief Flattening constructor given a contraction on the right hand
     * side and another expression on the left hand side.
     **/
    TensorContraction(
      IndexedTensor<F> *lhs, TensorContraction<F> *rhs
    ): factors(rhs.factors) {
      factors.push_back(lhs);
      // the factors from the rhs expression are now contained here
      rhs->factors.clear();
      delete rhs;
    }
    /**
     * \brief Constructor given two indexed tensors.
     **/
    TensorContraction(
      IndexedTensor<F> *lhs, IndexedTensor<F> *rhs
    ) {
      factors.push_back(lhs);
      factors.push_back(rhs);
    }
    /**
     * \brief Constructor given two general expressions.
     * This is currently not supported.
     **/
    TensorContraction(
      TensorExpression<F> *lhs, TensorExpression<F> *rhs
    ) {
      static_assert(
        StaticAssert<F>::False,
        "Only contractions of contractions or tensors supported."
      );
    }
    virtual ~TensorContraction() {
      // subexpressions are dependent entities: delete each factor
      for (auto factor(factors.begin()); factor != factors.end(); ++factor) {
        delete *factor;
      }
    }

    virtual void log() const {
      for (auto factor(factors.begin()); factor != factors.end(); ++factor) {
        (*factor)->log();
      }
      LOG(0, "TCC") << factors.size() << " tensors contracted" << std::endl;
    }

    std::vector<IndexedTensor<F> *> factors;
  };

  template <typename Lhs, typename Rhs>
  inline TensorContraction<typename Lhs::FieldType> &operator *(
    Lhs &A, Rhs &B
  ) {
    static_assert(
      TypeRelations<typename Lhs::FieldType, typename Rhs::FieldType>::Equals,
      "Only tensors of the same type can be contracted."
    );
    return *new TensorContraction<typename Lhs::FieldType>(&A, &B);
  }
}

#endif

