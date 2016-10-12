/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SCA_LAPACK_SINGULAR_VALUE_DECOMPOSITION_DEFINED
#define SCA_LAPACK_SINGULAR_VALUE_DECOMPOSITION_DEFINED

#include <math/Complex.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <extern/ScaLapack.hpp>

namespace cc4s {
  // base template
  template <typename F=double>
  class ScaLapackSingularValueDecomposition;

  // specialization for double
  template <>
  class ScaLapackSingularValueDecomposition<double> {
  public:
    ScaLapackSingularValueDecomposition(
      ScaLapackMatrix<double> const *A_,
      ScaLapackMatrix<double> *U_,
      ScaLapackMatrix<double> *VT_
    ):
      A(A_), U(U_), VT(VT_), workCount(-1), work(nullptr)
    {
      double optimalWork;
      int offset(1);
      int info;
      pdgesvd_(
        "V", "V",
        &A->lens[0], &A->lens[1],
        A->getLocalValues(), &offset, &offset, A->getDescriptor(),
        nullptr,
        U->getLocalValues(), &offset, &offset, U->getDescriptor(),
        VT->getLocalValues(), &offset, &offset, VT->getDescriptor(),
        &optimalWork, &workCount,
        &info
      );
      // TODO: check info
      workCount = static_cast<int>(optimalWork+0.5);
      // allocate work:
      work = new double[workCount];
    }

    ~ScaLapackSingularValueDecomposition() {
      if (work) delete[] work;
      work = nullptr;
    }

    void decompose(double *sigma) {
      int offset(1);
      int info;
      pdgesvd_(
        "V", "V",
        &A->lens[0], &A->lens[1],
        A->getLocalValues(), &offset, &offset, A->getDescriptor(),
        sigma,
        U->getLocalValues(), &offset, &offset, U->getDescriptor(),
        VT->getLocalValues(), &offset, &offset, VT->getDescriptor(),
        work, &workCount,
        &info
      );
      // TODO: check info
    }
  protected:
    ScaLapackMatrix<double> const *A;
    ScaLapackMatrix<double> *U, *VT;
    int workCount;
    double *work;
  };

  // specialization for complex
  template <>
  class ScaLapackSingularValueDecomposition<complex> {
  public:
    ScaLapackSingularValueDecomposition(
      ScaLapackMatrix<complex> const *A_,
      ScaLapackMatrix<complex> *U_,
      ScaLapackMatrix<complex> *VT_
    ):
      A(A_), U(U_), VT(VT_), workCount(-1), work(nullptr), realWork(nullptr)
    {
      complex optimalWork;
      double optimalRealWork;
      int offset(1);
      int info;
      pzgesvd_(
        "V", "V",
        &A->lens[0], &A->lens[1],
        A->getLocalValues(), &offset, &offset, A->getDescriptor(),
        nullptr,
        U->getLocalValues(), &offset, &offset, U->getDescriptor(),
        VT->getLocalValues(), &offset, &offset, VT->getDescriptor(),
        &optimalWork, &workCount, &optimalRealWork,
        &info
      );
      // TODO: check info
      workCount = static_cast<int>(std::real(optimalWork)+0.5);
      // allocate work:
      work = new complex[workCount];
      realWork = new double[static_cast<int64_t>(optimalRealWork+0.5)];
    }

    ~ScaLapackSingularValueDecomposition() {
      if (work) delete[] work;
      if (realWork) delete[] realWork;
      work = nullptr;
      realWork = nullptr;
    }
    void decompose(double *sigma) {
      int offset(1);
      int info;
      pzgesvd_(
        "V", "V",
        &A->lens[0], &A->lens[1],
        A->getLocalValues(), &offset, &offset, A->getDescriptor(),
        sigma,
        U->getLocalValues(), &offset, &offset, U->getDescriptor(),
        VT->getLocalValues(), &offset, &offset, VT->getDescriptor(),
        work, &workCount, realWork,
        &info
      );
      // TODO: check info
    }
  protected:
    ScaLapackMatrix<complex> const *A;
    ScaLapackMatrix<complex> *U, *VT;
    int workCount;
    complex *work;
    double *realWork;
  };
}

#endif

