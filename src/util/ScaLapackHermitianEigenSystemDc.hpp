/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SCA_LAPACK_HERMITIAN_EIGEN_SYSTEM_DC_DEFINED
#define SCA_LAPACK_HERMITIAN_EIGEN_STSTEM_DC_DEFINED

#include <math/Complex.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <extern/ScaLapack.hpp>
#include <util/Exception.hpp>

namespace cc4s {
  // base template
  template <typename F=double>
  class ScaLapackHermitianEigenSystemDc;

  // specialization for double
  template <>
  class ScaLapackHermitianEigenSystemDc<double> {
  public:
    ScaLapackHermitianEigenSystemDc(
      ScaLapackMatrix<double> const *A_,
      ScaLapackMatrix<double> *U_
    ):
      A(A_), U(U_),
      workCount(-1), integerWorkCount(-1),
      work(nullptr), integerWork(nullptr)
    {
      // TODO: check if matrix is quadratic
      double optimalWork;
      int optimalIntegerWork;
      int offset(1);
      int info;
      pdsyevd_(
        "V", "L",
        &A->lens[0],
        A->getLocalValues(), &offset, &offset, A->getDescriptor(),
        nullptr,
        U->getLocalValues(), &offset, &offset, U->getDescriptor(),
        &optimalWork, &workCount,
        &optimalIntegerWork, &integerWorkCount,
        &info
      );
      if (info != 0) {
        throw new EXCEPTION("ERROR: PDSYEVD initialization failed");
      }
      workCount = static_cast<int>(optimalWork+0.5);
      integerWorkCount = optimalIntegerWork;
      // allocate work:
      work = new double[workCount];
      integerWork = new int[integerWorkCount];
    }

    ~ScaLapackHermitianEigenSystemDc() {
      if (work) delete[] work;
      if (integerWork) delete[] integerWork;
      work = nullptr;
      integerWork = nullptr;
    }

    void solve(double *lambda) {
      int offset(1);
      int info;
      pdsyevd_(
        "V", "L",
        &A->lens[0],
        A->getLocalValues(), &offset, &offset, A->getDescriptor(),
        lambda,
        U->getLocalValues(), &offset, &offset, U->getDescriptor(),
        work, &workCount,
        integerWork, &integerWorkCount,
        &info
      );
      if (info != 0) {
        throw new EXCEPTION("ERROR: PDSYEVD diagonlization failed");
      }
    }
  protected:
    ScaLapackMatrix<double> const *A;
    ScaLapackMatrix<double> *U;
    int workCount, integerWorkCount;
    double *work;
    int *integerWork;
  };

  // specialization for complex
  template <>
  class ScaLapackHermitianEigenSystemDc<complex> {
  public:
    ScaLapackHermitianEigenSystemDc(
      ScaLapackMatrix<complex> const *A_,
      ScaLapackMatrix<complex> *U_
    ):
      A(A_), U(U_),
      workCount(-1), realWorkCount(-1), integerWorkCount(-1),
      work(nullptr), realWork(nullptr), integerWork(nullptr)
    {
      // TODO: check if matrix is quadratic
      complex optimalWork;
      double optimalRealWork;
      int optimalIntegerWork;
      int offset(1);
      int info;
      pzheevd_(
        "V", "L",
        &A->lens[0],
        A->getLocalValues(), &offset, &offset, A->getDescriptor(),
        nullptr,
        U->getLocalValues(), &offset, &offset, U->getDescriptor(),
        &optimalWork, &workCount,
        &optimalRealWork, &realWorkCount,
        &optimalIntegerWork, &integerWorkCount,
        &info
      );
      if (info != 0) {
        throw new EXCEPTION("ERROR: PZHEEVD initialization failed");
      }
      workCount = static_cast<int>(std::real(optimalWork)+0.5);
      realWorkCount = static_cast<int64_t>(optimalRealWork+0.5);
      integerWorkCount = optimalIntegerWork;
      // allocate work:
      work = new complex[workCount];
      // NOTE that realWorkCount actually refers to a number of complexes
      realWork = new double[2*realWorkCount];
      integerWork = new int[integerWorkCount];
    }

    ~ScaLapackHermitianEigenSystemDc() {
      if (work) delete[] work;
      if (realWork) delete[] realWork;
      if (integerWork) delete[] integerWork;
      work = nullptr;
      realWork = nullptr;
      integerWork = nullptr;
    }

    void solve(double *lambda) {
      int offset(1);
      int info;
      pzheevd_(
        "V", "L",
        &A->lens[0],
        A->getLocalValues(), &offset, &offset, A->getDescriptor(),
        lambda,
        U->getLocalValues(), &offset, &offset, U->getDescriptor(),
        work, &workCount,
        realWork, &realWorkCount,
        integerWork, &integerWorkCount,
        &info
      );
      if (info != 0) {
        throw new EXCEPTION("ERROR: PZHEEVD diagonlization failed");
      }
    }
  protected:
    ScaLapackMatrix<complex> const *A;
    ScaLapackMatrix<complex> *U;
    int workCount, realWorkCount, integerWorkCount;
    complex *work;
    double *realWork;
    int *integerWork;
  };
}

#endif

