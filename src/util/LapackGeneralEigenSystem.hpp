/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LAPACK_GENERAL_EIGEN_SYSTEM_DEFINED
#define LAPACK_GENERAL_EIGEN_SYSTEM_DEFINED

#include <math/Complex.hpp>
#include <util/LapackMatrix.hpp>
#include <util/Exception.hpp>
#include <extern/Lapack.hpp>
#include <util/Log.hpp>

#include <vector>

namespace cc4s {
  // base template
  template <typename F=double>
  class LapackGeneralEigenSystem;

/*
  // specialization for double
  template <>
  class LapackGeneralEigenSystem<double> {
  public:
    LapackGeneralEigenSystem(
      const LapackMatrix<double> &A_
    ):
      R(A_.getRows(), A_.getColumns()),
      L(A_.getRows(), A_.getColumns()),
      lambdas(A_.getRows())
    {
      if (A_.getRows() != A_.getColumns()) {
        throw EXCEPTION("EigenSystem requries a square matrix");
      }
      // copy A since it will be modified
      LapackMatrix<double> A(A_);
      int rows(A_.getRows());
      std::vector<double> lambdaReals(rows), lambdaImags(rows);
      double optimalWork;
      int workCount(-1);
      int info;
      dgeev_(
        "V", "V",
        &rows,
        A.getValues(), &rows,
        lambdaReals.data(), lambdaImags.data(),
        L.getValues(), &rows,
        R.getValues(), &rows,
        &optimalWork, &workCount,
        &info
      );
      // TODO: check info
      workCount = static_cast<int>(optimalWork+0.5);
      std::vector<double> work(workCount);
     
      dgeev_(
        "V", "V",
        &rows,
        A.getValues(), &rows,
        lambdaReals.data(), lambdaImags.data(),
        L.getValues(), &rows,
        R.getValues(), &rows,
        work.data(), &workCount,
        &info
      );
      // TODO: check info

      for (int i(0); i < rows; ++i) {
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
*/

  // specialization for complex
  template <>
  class LapackGeneralEigenSystem<complex> {
  public:
    LapackGeneralEigenSystem(
      const LapackMatrix<complex> &A_
    ):
      R(A_.getRows(), A_.getColumns()),
      L(A_.getRows(), A_.getColumns()),
      lambdas(A_.getRows())
    {
      if (A_.getRows() != A_.getColumns()) {
        throw EXCEPTION("EigenSystem requries a square matrix");
      }
      // copy A since it will be modified
      LapackMatrix<complex> A(A_);
      int rows(A_.getRows());
      std::vector<double> realWork(2*rows);
      complex optimalWork;
      int workCount(-1);
      int info;
      zgeev_(
        "V", "V",
        &rows,
        A.getValues(), &rows,
        lambdas.data(),
        L.getValues(), &rows,
        R.getValues(), &rows,
        &optimalWork, &workCount,
        realWork.data(),
        &info
      );
      // TODO: check info
      workCount = static_cast<int>(std::real(optimalWork)+0.5);
      std::vector<complex> work(workCount);
     
      zgeev_(
        "V", "V",
        &rows,
        A.getValues(), &rows,
        lambdas.data(),
        L.getValues(), &rows,
        R.getValues(), &rows,
        work.data(), &workCount,
        realWork.data(),
        &info
      );
      // TODO: check info
    }

    ~LapackGeneralEigenSystem() {
    }

    const std::vector<complex> &getEigenValues() const {
      return lambdas;
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

