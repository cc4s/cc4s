#ifndef SINGULAR_VALUE_DECOMPOSITION_DEFINED
#define SINGULAR_VALUE_DECOMPOSITION_DEFINED

#include <math/Complex.hpp>
#include <util/DryTensor.hpp>
#include <ctf.hpp>
#include <random>

namespace cc4s {
  template <typename F>
  class SingularValueDecomposition {
  public:
    SingularValueDecomposition(CTF::Matrix<F> const &matrix);
    CTF::Matrix<F> &get();

  protected:
    CTF::Matrix<F> inverse;
    void getProcessorGrid(
      int processors, int &processorRows, int &processorColumns
    );
  };

  template <typename F>
  class DrySingularValueDecomposition {
  public:
    DrySingularValueDecomposition(DryMatrix<F> const &matrix);
    DryMatrix<F> &get();

  protected:
    DryMatrix<F> inverse;
  };
}

#endif

