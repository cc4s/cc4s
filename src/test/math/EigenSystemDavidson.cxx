#include <test/Test.hpp>

#include <math/EigenSystemDavidson.hpp>
#include <util/LapackMatrix.hpp>
#include <util/LapackGeneralEigenSystem.hpp>
#include <math/Vector.hpp>

#include <util/Log.hpp>

#include <vector>
#include <random>
#include <utility>
#include <algorithm>

namespace cc4s {
  template <int D>
  class LinearMap {
  public:
    LinearMap(): A(D,D) {
    }

    Vector<complex,D> rightApply(const Vector<complex,D> &v) const {
      Vector<complex,D> w;
      for (int i(0); i < D; ++i) {
        for (int k(0); k < D; ++k) {
          w[i] += A(i,k) * v[k];
        }
      }
      return w;
    }

    LapackMatrix<complex> A;
  };

  template <int D>
  class DiagonalDominantPreconditioner {
  public:
    DiagonalDominantPreconditioner(const LinearMap<D> &f) {
      for (int i(0); i < D; ++i) {
        diagonal[i] = f.A(i,i);
      }
    }

    std::vector<Vector<complex,D>> getInitialBasis(const int N) const {
      // return unit vectors to lowest elemens in diagonal
      Vector<complex,D> d(diagonal);
      std::vector<Vector<complex,D>> basis(N);
      double maxReal(std::real(d[0]));
      for (int i(1); i < D; ++i) maxReal = std::max( std::real(d[i]), maxReal );
      for (int b(0); b < N; ++b) {
        double minReal(std::real(d[0]));
        int minIndex(0);
        for (int i(1); i < D; ++i) {
          if (std::real(d[i]) < minReal) {
            minReal = std::real(d[i]);
            minIndex = i;
          }
        }
        // set to maximum value not to encounter it again
        d[minIndex] = maxReal;
        Vector<complex,D> v;
        v[minIndex] = 1.0;
        basis[b] = v;
      }
      return basis;
    }

    Vector<complex,D> getCorrection(
      const complex lambda, const Vector<complex,D> &residuum
    ) const {
      Vector<complex,D> w;
      // compute ((lambda * id - Diag(diagonal))^-1) . residuum
      for (int k(0); k < D; ++k) {
        w[k] = residuum[k] / (lambda - diagonal[k]);
      }
      return w;
    }

    Vector<complex,D> diagonal;
  };
}

using namespace cc4s;

TEST_CASE( "EigenSystemDavidson", "[math]" ) {
  constexpr int N(256);
  std::mt19937 random;
  std::normal_distribution<double> normalDistribution(0.0, 1.0);
  LinearMap<N> f;
  for (int i(0); i < N; ++i) {
    for (int j(0); j < N; ++j) {
      f.A(i,j).real(normalDistribution(random));
      f.A(i,j).imag(0);
      //f.A(i,j).imag(normalDistribution(random));
    }
    // make it diagonally dominant
    f.A(i,i) -= 100.0;
  }
  f.A(10,10) = -100.0;

  // Find the eigenvalues using the Lapack eigen solver
  LapackGeneralEigenSystem<complex> eigenSystem(f.A);

  EigenSystemDavidson<Vector<complex,N>> eigenSystemDavidson(
    f, 4, DiagonalDominantPreconditioner<N>(f), 1E-8, 64
  );

  for (unsigned int k(0); k < eigenSystemDavidson.getEigenValues().size(); ++k) {
    std::stringstream rowStream;
    for (int i(0); i < eigenSystem.getRightEigenVectors().getRows(); ++i) {
      rowStream << "\t" << eigenSystem.getRightEigenVectors()(i,k);
    }
    LOG(0, "LapackGeneralEigenSystem") <<
      "EigenValue[" << k << "]=" <<
      eigenSystem.getEigenValues()[k] << std::endl;
    LOG(0, "LapackGeneralEigenSystem") <<
      "EigenValueDavison[" << k << "]=" <<
      eigenSystemDavidson.getEigenValues()[k] << std::endl;
/*
    LOG(0, "LapackGeneralEigenSystem") <<
      "EigenVector[" << k << "]=" << rowStream.str() << std::endl;
    LOG(0, "LapackGeneralEigenSystem") <<
      "EigenVectorDavidson[" << k << "]=" <<
      eigenSystemDavidson.getRightEigenVectors()[k] << std::endl;
*/
  };
}
