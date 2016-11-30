#ifndef PSEUDO_INVERSE_SVD_DEFINED
#define PSEUDO_INVERSE_SVD_DEFINED

#include <math/Complex.hpp>
#include <tcc/DryTensor.hpp>

#include <ctf.hpp>
#include <random>

namespace cc4s {
  template <typename F>
  class PseudoInverseSvd {
  public:
    PseudoInverseSvd(CTF::Matrix<F> &A, double epsilon = 1e-12);
    CTF::Matrix<F> &get();

  protected:
    CTF::Matrix<F> inverse;
  };

  template <typename F>
  class DryPseudoInverseSvd {
  public:
    DryPseudoInverseSvd(DryMatrix<F> const &matrix);
    DryMatrix<F> &get();

  protected:
    DryMatrix<F> inverse;
  };
}

#endif

