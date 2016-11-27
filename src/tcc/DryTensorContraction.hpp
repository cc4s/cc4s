/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_TENSOR_CONTRACTION_DEFINED
#define DRY_TENSOR_CONTRACTION_DEFINED

#include <tcc/DryTensorExpression.hpp>
#include <util/Log.hpp>
#include <vector>

namespace cc4s {
  template <typename F>
  class DryTensorContraction: public DryTensorExpression<F> {
  public:
    /**
     * \brief Flattening constructor given two contractions.
     **/
    DryTensorContraction(
      DryTensorContraction<F> *lhs, DryTensorContraction<F> *rhs
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
    DryTensorContraction(
      DryTensorContraction<F> *lhs, DryTensorExpression<F> *rhs
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
    DryTensorContraction(
      DryTensorExpression<F> *lhs, DryTensorContraction<F> *rhs
    ): factors(rhs.factors) {
      factors.push_back(lhs);
      // the factors from the rhs expression are now contained here
      rhs->factors.clear();
      delete rhs;
    }
    /**
     * \brief Constructor given two general expressions.
     **/
    DryTensorContraction(
      DryTensorExpression<F> *lhs, DryTensorExpression<F> *rhs
    ) {
      factors.push_back(lhs);
      factors.push_back(rhs);
    }
    virtual ~DryTensorContraction() {
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

  protected:
    std::vector<DryTensorExpression<F> *> factors;
  };

  template <typename F>
  DryTensorContraction<F> &operator *(
    DryTensorContraction<F> &A, DryTensorContraction<F> &B
  ) {
    return *new DryTensorContraction<F>(&A, &B);
  }
  template <typename F>
  DryTensorContraction<F> &operator *(
    DryTensorContraction<F> &A, DryTensorExpression<F> &B
  ) {
    return *new DryTensorContraction<F>(&A, &B);
  }
  template <typename F>
  DryTensorContraction<F> &operator *(
    DryTensorExpression<F> &A, DryTensorContraction<F> &B
  ) {
    return *new DryTensorContraction<F>(&A, &B);
  }
  template <typename F>
  DryTensorContraction<F> &operator *(
    DryTensorExpression<F> &A, DryTensorExpression<F> &B
  ) {
    return *new DryTensorContraction<F>(&A, &B);
  }
}

#endif

