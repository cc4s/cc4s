/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRY_TENSOR_CONTRACTION_DEFINED
#define DRY_TENSOR_CONTRACTION_DEFINED

#include <util/DryTensorTerm.hpp>
#include <util/Log.hpp>

namespace cc4s {
  template <typename F=double>
  class DryTensorContraction: public cc4s::DryTensorTerm<F> {
  public:
    DryTensorContraction(
      cc4s::DryTensorTerm<F> const &A_, cc4s::DryTensorTerm<F> const &B_
    ): A(&A_), B(&B_) {
    }
    virtual ~DryTensorContraction() {
    }

    virtual void log() const {
      A->log();
      B->log();
      LOG(0, "TCC") << "contracted" << std::endl;
    }

  protected:
    cc4s::DryTensorTerm<F> const *A, *B;
  };

  template <typename F>
  cc4s::DryTensorContraction<F> operator *(
    cc4s::DryTensorTerm<F> const &A, cc4s::DryTensorTerm<F> const &B
  ) {
    return cc4s::DryTensorContraction<F>(A, B);
  }
}

#endif

