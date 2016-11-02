/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SCA_LAPACK_HERMITIAN_EIGEN_SYSTEM_DEFINED
#define SCA_LAPACK_HERMITIAN_EIGEN_STSTEM_DEFINED

#include <math/Complex.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <extern/ScaLapack.hpp>

namespace cc4s {
  // base template
  template <typename F=double>
  class ScaLapackHermitianEigenSystem;

  // specialization for double
  template <>
  class ScaLapackHermitianEigenSystem<double> {
  public:
    ScaLapackHermitianEigenSystem(
      ScaLapackMatrix<double> const *A_,
      ScaLapackMatrix<double> *U_
    ):
      A(A_), U(U_), workCount(-1), work(nullptr)
    {
      // TODO: check if matrix is quadratic
      double optimalWork;
      int offset(1);
      int info;
      pdsyev_(
        "V", "L",
        &A->lens[0],
        A->getLocalValues(), &offset, &offset, A->getDescriptor(),
        nullptr,
        U->getLocalValues(), &offset, &offset, U->getDescriptor(),
        &optimalWork, &workCount,
        &info
      );
      // TODO: check info
      workCount = static_cast<int>(optimalWork+0.5);
      // allocate work:
      work = new double[workCount];
    }

    ~ScaLapackHermitianEigenSystem() {
      if (work) delete[] work;
      work = nullptr;
    }

    void solve(double *lambda) {
      int offset(1);
      int info;
      pdsyev_(
        "V", "L",
        &A->lens[0],
        A->getLocalValues(), &offset, &offset, A->getDescriptor(),
        lambda,
        U->getLocalValues(), &offset, &offset, U->getDescriptor(),
        work, &workCount,
        &info
      );
      // TODO: check info
    }
  protected:
    ScaLapackMatrix<double> const *A;
    ScaLapackMatrix<double> *U;
    int workCount;
    double *work;
  };

  // specialization for complex
  template <>
  class ScaLapackHermitianEigenSystem<complex> {
  public:
    ScaLapackHermitianEigenSystem(
      ScaLapackMatrix<complex> const *A_,
      ScaLapackMatrix<complex> *U_
    ):
      A(A_), U(U_),
      workCount(-1), realWorkCount(-1), work(nullptr), realWork(nullptr)
    {
      // TODO: check if matrix is quadratic
      complex optimalWork;
      double optimalRealWork;
      int offset(1);
      int info;
      pzheev_(
        "V", "L",
        &A->lens[0],
        A->getLocalValues(), &offset, &offset, A->getDescriptor(),
        nullptr,
        U->getLocalValues(), &offset, &offset, U->getDescriptor(),
        &optimalWork, &workCount, &optimalRealWork, &realWorkCount,
        &info
      );
      // TODO: check info
      workCount = static_cast<int>(std::real(optimalWork)+0.5);
      realWorkCount = static_cast<int64_t>(optimalRealWork+0.5);
      // allocate work:
      work = new complex[workCount];
      // NOTE that realWorkCount actually refers to a number of complexes
      realWork = new double[2*realWorkCount];
    }

    ~ScaLapackHermitianEigenSystem() {
      if (work) delete[] work;
      if (realWork) delete[] realWork;
      work = nullptr;
      realWork = nullptr;
    }

    void solve(double *lambda) {
      int offset(1);
      int info;
      pzheev_(
        "V", "L",
        &A->lens[0],
        A->getLocalValues(), &offset, &offset, A->getDescriptor(),
        lambda,
        U->getLocalValues(), &offset, &offset, U->getDescriptor(),
        work, &workCount, realWork, &realWorkCount,
        &info
      );
      // TODO: check info
    }
  protected:
    ScaLapackMatrix<complex> const *A;
    ScaLapackMatrix<complex> *U;
    int workCount, realWorkCount;
    complex *work;
    double *realWork;
  };
}

#endif

