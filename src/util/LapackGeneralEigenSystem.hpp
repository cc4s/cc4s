/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LAPACK_GENERAL_EIGEN_SYSTEM_DEFINED
#define LAPACK_GENERAL_EIGEN_SYSTEM_DEFINED

#include <math/Complex.hpp>
#include <util/LapackMatrix.hpp>

namespace cc4s {
  // base template
  template <typename F=double>
  class LapackGeneralEigenSystem;

  // specialization for double
  template <>
  class LapackGeneralEigenSystem<double> {
  public:
    LapackGeneralEigenSystem(
      const LapackMatrix<double> &A_,
      LapackMatrix<double> *R_
    ):
      A(A_), R(R_), workCount(-1), work(nullptr)
    {
      // TODO: check if matrix is quadratic
      double optimalWork;
      int info;
      int one(1);
      dgeev_(
        "N", "V", // compute only right eigenvectors
        &A.rows,
        A.getValues(), &A.rows,
        nullptr, nullptr,
        nullptr, &one,
        &R->getValues(), &R->rows,
        &optimalWork, &workCount,
        &info
      );
      // TODO: check info
      workCount = static_cast<int>(optimalWork+0.5);
      // allocate work:
      work = new double[workCount];
    }

    ~LapackGeneralEigenSystem() {
      if (work) delete[] work;
      work = nullptr;
    }

    void solve(complex *lambda) {
      int offset(1);
      int info;
      pdsyev_(
        "V", "L",
        &A->lens[0],
        A->getLocalValues(), &offset, &offset, A->getDescriptor(),
        lambda,
        R->getLocalValues(), &offset, &offset, U->getDescriptor(),
        work, &workCount,
        &info
      );
      // TODO: check info
    }
  protected:
    LapackMatrix<double> A;
    LapackMatrix<double> *R;
    int workCount;
    double *work;
  };

  // TODO: write complex version
  // specialization for complex
}

#endif

