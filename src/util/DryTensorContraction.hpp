/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_TENSOR_CONTRACTION_DEFINED
#define DRY_TENSOR_CONTRACTION_DEFINED

#include <util/DryTensorExpression.hpp>
#include <util/Log.hpp>
#include <vector>

namespace cc4s {
  template <typename F=double>
  class DryTensorContraction: public cc4s::DryTensorExpression<F> {
  public:
    /**
     * \brief Flattening constructor given two contractions.
     **/
    DryTensorContraction(
      DryTensorContraction<F> const &lhs, DryTensorContraction<F> const &rhs
    ): factors(lhs.factors) {
      factors.insert(factors.end(), rhs.factors.begin(), rhs.factors.end());
    }
    /**
     * \brief Flattening constructor given a contraction on the left hand
     * side and another expression on the right hand side.
     **/
    DryTensorContraction(
      DryTensorContraction<F> const &lhs, DryTensorExpression<F> const &rhs
    ): factors(lhs.factors) {
      factors.push_back(&rhs);
    }
    /**
     * \brief Flattening constructor given a contraction on the right hand
     * side and another expression on the left hand side.
     **/
    DryTensorContraction(
      DryTensorExpression<F> const &lhs, DryTensorContraction<F> const &rhs
    ): factors(rhs.factors) {
      factors.push_back(&lhs);
    }
    /**
     * \brief Constructor given two general expressions.
     **/
    DryTensorContraction(
      DryTensorExpression<F> const &lhs, DryTensorExpression<F> const &rhs
    ) {
      factors.push_back(&lhs);
      factors.push_back(&rhs);
    }
    virtual ~DryTensorContraction() {
    }

    virtual void log() const {
      for (auto factor(factors.begin()); factor != factors.end(); ++factor) {
        (*factor)->log();
      }
      LOG(0, "TCC") << factors.size() << " tensors contracted" << std::endl;
    }

  protected:
    std::vector<cc4s::DryTensorExpression<F> const *> factors;
  };

  template <typename F>
  cc4s::DryTensorContraction<F> operator *(
    cc4s::DryTensorContraction<F> const &A, cc4s::DryTensorContraction<F> const &B
  ) {
    return cc4s::DryTensorContraction<F>(A, B);
  }
  template <typename F>
  cc4s::DryTensorContraction<F> operator *(
    cc4s::DryTensorContraction<F> const &A, cc4s::DryTensorExpression<F> const &B
  ) {
    return cc4s::DryTensorContraction<F>(A, B);
  }
  template <typename F>
  cc4s::DryTensorContraction<F> operator *(
    cc4s::DryTensorExpression<F> const &A, cc4s::DryTensorContraction<F> const &B
  ) {
    return cc4s::DryTensorContraction<F>(A, B);
  }
  template <typename F>
  cc4s::DryTensorContraction<F> operator *(
    cc4s::DryTensorExpression<F> const &A, cc4s::DryTensorExpression<F> const &B
  ) {
    return cc4s::DryTensorContraction<F>(A, B);
  }
}

#endif

