/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LAPACK_GENERAL_EIGEN_SYSTEM_DEFINED
#define LAPACK_GENERAL_EIGEN_SYSTEM_DEFINED

#include <math/Complex.hpp>
#include <util/LapackMatrix.hpp>

#include <vector>

namespace cc4s {
  // base template
  template <typename F=double>
  class LapackGeneralEigenSystem;

  // specialization for double
  template <>
  class LapackGeneralEigenSystem<double> {
  public:
    LapackGeneralEigenSystem(
      const LapackMatrix<double> &A_
    ):
      R(),
      lambdas(A_.getRows()),
    {
      if (A_.getRows() != A_.getColumns()) {
        throw EXCEPTION("EigenSystem requries a square matrix");
      }
      // copy A since it will be modified
      LapackMatrix<double> A(A_);
      std::vector<double> lambdaReals(A_.getRows()), lambdaImags(A_.getRows());
      double optimalWork;
      int workCount(-1);
      int info;
      dgeev_(
        "N", "V", // compute only right eigenvectors
        &A.getRows(),
        A.data(), &A.getRows(),
        lambdaReals.data(), lambdaImags.data(),
        nullptr, &one,
        &R->getValues(), &R->getRows(),
        &optimalWork, &workCount,
        &info
      );
      // TODO: check info
      workCount = static_cast<int>(optimalWork+0.5);
      std::vector<double> work(workCount);
     
      dgeev_(
        "N", "V", // compute only right eigenvectors
        &A.getRows(),
        A.data(), &A.getRows(),
        lambdaReals.data(), lambdaImags.data(),
        nullptr, &one,
        &R->getValues(), &R->getRows(),
        &work.data(), &workCount,
        &info
      );
      // TODO: check info

      for (int i(0); i < A.getRows(); ++i) {
        lambdas[i] = complex(lambdaReals[i], lambdaImags[i]);
        if (std::abs(lambdaImags[i]) > 1e-8*std::abs(lambdaReals[i])) {
          // TODO: decode eigenvectors to complex eigenvalues
        }
      }
    }

    ~LapackGeneralEigenSystem() {
    }

    const std::vector<complex> &getEigenValues() const {
      return lambda;
    }

    const LapackMatrix<complex> &getRightEigenVectors() const {
      return R;
    }

    const LapackMatrix<complex> &getLeftEigenVectors() const {
      return L;
    }

  protected:
    LapackMatrix<complex> R, L;
    std::vector<complex> lambdas;
  };

  // specialization for complex
  template <>
  class LapackGeneralEigenSystem<complex> {
  public:
    LapackGeneralEigenSystem(
      const LapackMatrix<complex> &A_
    ):
      R(),
      lambdas(A_.getRows()),
    {
      if (A_.getRows() != A_.getColumns()) {
        throw EXCEPTION("EigenSystem requries a square matrix");
      }
      // copy A since it will be modified
      LapackMatrix<complex> A(A_);
      std::vector<complex> lambdaReals(A_.getRows()), lambdaImags(A_.getRows());
      complex optimalWork;
      int workCount(-1);
      int info;
      zgeev_(
        "N", "V", // compute only right eigenvectors
        &A.getRows(),
        A.data(), &A.getRows(),
        lambdas.data(),
        nullptr, &one,
        &R->getValues(), &R->getRows(),
        &optimalWork, &workCount,
        &info
      );
      // TODO: check info
      workCount = static_cast<int>(std::real(optimalWork)+0.5);
      std::vector<complex> work(workCount);
      std::vector<double> realWork(2*A.getRows());
     
      zgeev_(
        "N", "V", // compute only right eigenvectors
        &A.getRows(),
        A.data(), &A.getRows(),
        lambdas.data(),
        nullptr, &one,
        &R->getValues(), &R->getRows(),
        &work.data(), &workCount,
        &realWorl.data(),
        &info
      );
      // TODO: check info
    }

    ~LapackGeneralEigenSystem() {
    }

    const std::vector<complex> &getEigenValues() const {
      return lambda;
    }

    const LapackMatrix<complex> &getRightEigenVectors() const {
      return R;
    }

    const LapackMatrix<complex> &getLeftEigenVectors() const {
      return L;
    }

  protected:
    LapackMatrix<complex> R, L;
    std::vector<complex> lambdas;
  };
}

#endif

