/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
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


  template <>
  class LapackGeneralEigenSystem<double> {
  public:
    class EigenValueComparator {
    public:
      bool operator ()(
        const std::pair<int, double> &a,
        const std::pair<int, double> &b
      ) {
        double diff(b.second-a.second);
        double magnitude( std::abs(a.second)+std::abs(b.second) );
        if (std::real(diff) > +1E-13*magnitude) return true;
        if (std::real(diff) < -1E-13*magnitude) return false;
        return a.first < b.first;
      }
    };
  };

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

      // sort eigenvalues
      std::vector<std::pair<int, complex>> sortedEigenValues(rows);
      for (int i(0); i < rows; ++i) {
        sortedEigenValues[i] = std::pair<int, complex>(i, lambdas[i]);
      }
      std::sort(
        sortedEigenValues.begin(), sortedEigenValues.end(),
        EigenValueComparator()
      );
      // order eigenvectors and returned eigenvalues in the same way
      orderEigenVectors(sortedEigenValues, R);
      orderEigenVectors(sortedEigenValues, L);
      for (int i(0); i < rows; ++i) {
        lambdas[i] = sortedEigenValues[i].second;
      }
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

    double rightEigenError(const LapackMatrix<complex> &A) {
      double error(0);
      for (int i(0); i < A.getRows(); ++i) {
        for (int k(0); k < A.getRows(); ++k) {
          complex element(0);
          for (int j(0); j < A.getRows(); ++j) {
            element += A(i,j) * R(j,k);
          }
          element -= lambdas[k] * R(i,k);
          error += std::real(element*std::conj(element));
        }
      }
      return error;
    }
    double leftEigenError(const LapackMatrix<complex> &A) {
      double error(0);
      for (int j(0); j < A.getRows(); ++j) {
        for (int k(0); k < A.getRows(); ++k) {
          complex element(0);
          for (int i(0); i < A.getRows(); ++i) {
            element += std::conj(L(i,k)) * A(i,j);
          }
          element -= std::conj(L(j,k)) * lambdas[k];
          error += std::real(element*std::conj(element));
        }
      }
      return error;
    }

    double biorthogonalError() {
      double error(0);
      for (int i(0); i < R.getRows(); ++i) {
        for (int j(0); j < R.getRows(); ++j) {
          complex element(0);
          for (int k(0); k < R.getRows(); ++k) {
            element += std::conj(R(i,k)) * R(j,k);
          }
          element -= i == j ? 1.0 : 0.0;
          error += i == j ? 0.0 : std::real(element*std::conj(element));
        }
      }
      return error;
    }

    class EigenValueComparator {
    public:
      bool operator ()(
        const std::pair<int, complex> &a,
        const std::pair<int, complex> &b
      ) {
        // FIXME: Move the inf code into particular implementations
        double inf = std::numeric_limits<double>::infinity();
        complex A(
          std::abs(a.second) < 1E-4 ? *(new complex(inf,inf)) : a.second
        );
        complex B(
          std::abs(b.second) < 1E-4 ? *(new complex(inf,inf)) : b.second
        );
        complex diff(B-A);
        double magnitude( std::abs(a.second)+std::abs(b.second) );
        if (std::real(diff) > +1E-13*magnitude) return true;
        if (std::real(diff) < -1E-13*magnitude) return false;
        if (std::imag(diff) > +1E-13*magnitude) return false;
        if (std::imag(diff) < -1E-13*magnitude) return true;
        return a.first < b.first;
      }
    };

  protected:
    LapackMatrix<complex> R, L;
    std::vector<complex> lambdas;

    void orderEigenVectors(
      const std::vector<std::pair<int, complex>> &sortedEigenValues,
      LapackMatrix<complex> &U
    ) {
      LapackMatrix<complex> unsortedU(U);
      for (int j(0); j < U.getRows(); ++j) {
        int unsortedJ( sortedEigenValues[j].first );
        for (int i(0); i < U.getRows(); ++i) {
          U(i,j) = unsortedU(i,unsortedJ);
        }
      }
    }
  };
}

#endif

