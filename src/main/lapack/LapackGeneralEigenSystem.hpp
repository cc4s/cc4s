/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef LAPACK_GENERAL_EIGEN_SYSTEM_DEFINED
#define LAPACK_GENERAL_EIGEN_SYSTEM_DEFINED

#include <extern/Lapack.hpp>
#include <Real.hpp>
#include <Complex.hpp>
#include <LapackMatrix.hpp>
#include <Exception.hpp>
#include <Log.hpp>

#include <vector>

namespace cc4s {
  // base template
  template <typename F=Real<>>
  class LapackGeneralEigenSystem;


  template <>
  class LapackGeneralEigenSystem<Real<64>> {
  public:
    class EigenValueComparator {
    public:
      bool operator ()(
        const std::pair<int, Real<64>> &a,
        const std::pair<int, Real<64>> &b
      ) {
        Real<64> diff(b.second-a.second);
        Real<64> magnitude( abs(a.second)+abs(b.second) );
        if (real(diff) > +1E-13*magnitude) return true;
        if (real(diff) < -1E-13*magnitude) return false;
        return a.first < b.first;
      }
    };
  };

/*
  // specialization for Real<64>
  template <>
  class LapackGeneralEigenSystem<Real<64>> {
  public:
    LapackGeneralEigenSystem(
      const LapackMatrix<Real<64>> &A_
    ):
      R(A_.getRows(), A_.getColumns()),
      L(A_.getRows(), A_.getColumns()),
      lambdas(A_.getRows())
    {
      if (A_.getRows() != A_.getColumns()) {
        throw EXCEPTION("EigenSystem requries a square matrix");
      }
      // copy A since it will be modified
      LapackMatrix<Real<64>> A(A_);
      int rows(A_.getRows());
      std::vector<Real<64>> lambdaReals(rows), lambdaImags(rows);
      Real<64> optimalWork;
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
      std::vector<Real<64>> work(workCount);
     
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
        lambdas[i] = Complex<64>(lambdaReals[i], lambdaImags[i]);
        if (abs(lambdaImags[i]) > 1e-8*abs(lambdaReals[i])) {
          // TODO: decode eigenvectors to complex eigenvalues
        }
      }
    }

    ~LapackGeneralEigenSystem() {
    }

    const std::vector<Complex<64>> &getEigenValues() const {
      return lambda;
    }

    const LapackMatrix<Complex<64>> &getRightEigenVectors() const {
      return R;
    }

    const LapackMatrix<Complex<64>> &getLeftEigenVectors() const {
      return L;
    }

  protected:
    LapackMatrix<Complex<64>> R, L;
    std::vector<Complex<64>> lambdas;
  };
*/

  // specialization for complex
  template <>
  class LapackGeneralEigenSystem<Complex<64>> {
  public:
    LapackGeneralEigenSystem(
      const LapackMatrix<Complex<64>> &A_,
      bool computeRightEigenvectors = true,
      bool computeLeftEigenvectors = false
    ):
      R(
        computeRightEigenvectors ?
          NEW(LapackMatrix<Complex<64>>, A_.getRows(), A_.getColumns()) : nullptr
      ),
      L(
        computeLeftEigenvectors ?
          NEW(LapackMatrix<Complex<64>>, A_.getRows(), A_.getColumns()) : nullptr
      ),
      lambdas(A_.getRows())
    {
      if (A_.getRows() != A_.getColumns()) {
        throw EXCEPTION("EigenSystem requries a square matrix");
      }
      // copy A since it will be modified
      LapackMatrix<Complex<64>> A(A_);
      int rows(A_.getRows());
      std::vector<Real<64>> realWork(2*rows);
      Complex<64> optimalWork;
      int workCount(-1);
      int info;
      zgeev_(
        computeLeftEigenvectors ? "V" : "N",
        computeRightEigenvectors ? "V" : "N",
        &rows,
        A.getValues(), &rows,
        lambdas.data(),
        computeLeftEigenvectors ? L->getValues() : nullptr, &rows,
        computeRightEigenvectors ? R->getValues() : nullptr, &rows,
        &optimalWork, &workCount,
        realWork.data(),
        &info
      );
      // TODO: check info
      workCount = static_cast<int>(real(optimalWork)+0.5);
      std::vector<Complex<64>> work(workCount);

      zgeev_(
        computeLeftEigenvectors ? "V" : "N",
        computeRightEigenvectors ? "V" : "N",
        &rows,
        A.getValues(), &rows,
        lambdas.data(),
        computeLeftEigenvectors ? L->getValues() : nullptr, &rows,
        computeRightEigenvectors ? R->getValues() : nullptr, &rows,
        work.data(), &workCount,
        realWork.data(),
        &info
      );
      // TODO: check info

      // sort eigenvalues
      std::vector<std::pair<int, Complex<64>>> sortedEigenValues(rows);
      for (int i(0); i < rows; ++i) {
        sortedEigenValues[i] = std::pair<int, Complex<64>>(i, lambdas[i]);
      }
      std::sort(
        sortedEigenValues.begin(), sortedEigenValues.end(),
        EigenValueComparator()
      );
      // order eigenvectors and returned eigenvalues in the same way
      if (computeRightEigenvectors) {
        orderEigenVectors(sortedEigenValues, *R);
      }
      if (computeLeftEigenvectors) {
        orderEigenVectors(sortedEigenValues, *L);
      }
      for (int i(0); i < rows; ++i) {
        lambdas[i] = sortedEigenValues[i].second;
      }
    }

    ~LapackGeneralEigenSystem() {
    }

    const std::vector<Complex<64>> &getEigenValues() const {
      return lambdas;
    }

    const LapackMatrix<Complex<64>> &getRightEigenVectors() const {
      if (!R) {
        throw new EXCEPTION("Right eigenvectors were not computed in constructor.");
      }
      return *R;
    }

    const LapackMatrix<Complex<64>> &getLeftEigenVectors() const {
      if (!L) {
        throw new EXCEPTION("Left eigenvectors were not computed in constructor.");
      }
      return *L;
    }

    Real<64> rightEigenError(const LapackMatrix<Complex<64>> &A) {
      Real<64> error(0);
      for (int i(0); i < A.getRows(); ++i) {
        for (int k(0); k < A.getRows(); ++k) {
          Complex<64> element(0);
          for (int j(0); j < A.getRows(); ++j) {
            element += A(i,j) * (*R)(j,k);
          }
          element -= lambdas[k] * (*R)(i,k);
          error += real(element*conj(element));
        }
      }
      return error;
    }
    Real<64> leftEigenError(const LapackMatrix<Complex<64>> &A) {
      Real<64> error(0);
      for (int j(0); j < A.getRows(); ++j) {
        for (int k(0); k < A.getRows(); ++k) {
          Complex<64> element(0);
          for (int i(0); i < A.getRows(); ++i) {
            element += conj((*L)(i,k)) * A(i,j);
          }
          element -= conj((*L)(j,k)) * lambdas[k];
          error += real(element*conj(element));
        }
      }
      return error;
    }

    Real<64> biorthogonalError() {
      Real<64> error(0);
      for (int i(0); i < R->getRows(); ++i) {
        for (int j(0); j < R->getRows(); ++j) {
          Complex<64> element(0);
          for (int k(0); k < R->getRows(); ++k) {
            element += conj((*R)(i,k)) * (*R)(j,k);
          }
          element -= i == j ? 1.0 : 0.0;
          error += i == j ? 0.0 : real(element*conj(element));
        }
      }
      return error;
    }

    class EigenValueComparator {
    public:
      bool operator ()(
        const std::pair<int, Complex<64>> &a,
        const std::pair<int, Complex<64>> &b
      ) {
        // FIXME: Move the inf code into particular implementations
        Real<64> inf = std::numeric_limits<Real<64>>::infinity();
        Complex<64> A(
          abs(a.second) < 1E-4 ? *(new Complex<64>(inf,inf)) : a.second
        );
        Complex<64> B(
          abs(b.second) < 1E-4 ? *(new Complex<64>(inf,inf)) : b.second
        );
        Complex<64> diff(B-A);
        Real<64> magnitude( abs(a.second)+abs(b.second) );
        if (real(diff) > +1E-13*magnitude) return true;
        if (real(diff) < -1E-13*magnitude) return false;
        if (imag(diff) > +1E-13*magnitude) return false;
        if (imag(diff) < -1E-13*magnitude) return true;
        return a.first < b.first;
      }
    };

  protected:
    Ptr<LapackMatrix<Complex<64>>> R, L;
    std::vector<Complex<64>> lambdas;

    void orderEigenVectors(
      const std::vector<std::pair<int, Complex<64>>> &sortedEigenValues,
      LapackMatrix<Complex<64>> &U
    ) {
      LapackMatrix<Complex<64>> unsortedU(U);
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

