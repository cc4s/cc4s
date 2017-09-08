#ifndef EIGEN_SYSTEM_DAVIDSON_DEFINED
#define EIGEN_SYSTEM_DAVIDSON_DEFINED

#include <util/LapackMatrix.hpp>
#include <util/LapackGeneralEigenSystem.hpp>

#include <vector>
#include <utility>

namespace cc4s {
  template <typename V>
  class EigenSystemDavidson {
  public:
    typedef typename V::FieldType F;

    /**
     * \brief ...
     * \param[in] h object representing the matrix whose eigen system is sought
     * offering the following method:
     * V rightApply(V &v);
     * returning the action of h on v standing right of h.
     * \param[in] eigenVectorsCount the number of eigenvalues and vectors to be
     * computed.
     * \param[in] p object representing the preconditioner
     * offering the following method:
     * std::vector<V> getInitialBasis(int eigenVectorsCount);
     * returning an initial guess the basis consisting of eigenVectorsCount
     * vectors.
     * V getCorrection(const complex eigenValue, V &residuum);
     * computing the estimated correction for the k-th eigenvector given its
     * eigenvalue and the residuum of the current k-th estimated eigenvector.
     * \param[in] tolerance the targeted relative tolerance to be met
     * by all residua.
     * \param[in] maxBasisSize the maximum allowed number of vectors
     * representing eigenvectors.
     **/
    template <typename H, typename P>
    EigenSystemDavidson(
      H &h,
      const int eigenVectorsCount,
      P &p,
      const double tolerance = 1E-14,
      const unsigned int maxBasisSize = 1000
    ):
      eigenValues(eigenVectorsCount),
      rightEigenVectors( p.getInitialBasis(eigenVectorsCount) )
    {
      // get inital estimates for rEV = initial B matrix
      std::vector<V> basis( rightEigenVectors );

      // begin convergence loop
      double rms;
      int iterationCount(0);
      do {
        LOG(1,"Davidson") << "iteration=" << (iterationCount+1) << std::endl;
        // compute reduced H by projection onto subspace spanned by basis
        LapackMatrix<complex> reducedH(basis.size(), basis.size());
        for (unsigned int j(0); j < basis.size(); ++j) {
          V HBj( h.rightApply(basis[j]) );
          for (unsigned int i(0); i < basis.size(); ++i) {
            reducedH(i,j) = basis[i].dot(HBj);
          }
        }

        // compute K lowest reduced eigenvalues and vectors of reduced H
        LapackMatrix<complex> reducedEigenVectors(basis.size(), basis.size());
        LapackGeneralEigenSystem<complex> reducedEigenSystem(reducedH);

        // begin basis extension loop for each k
        rms = 0.0;
        for (unsigned int k(0); k < eigenValues.size(); ++k) {
          // get estimated eigenvalue
          eigenValues[k] = reducedEigenSystem.getEigenValues()[k];

          // compute estimated eigenvector by expansion in basis
          rightEigenVectors[k] *= F(0);
          for (int b(0); b < reducedH.getColumns(); ++b) {
            V scaledBase(
              basis[b] * ComplexTraits<F>::convert(
                reducedEigenSystem.getRightEigenVectors()(b,k)
              )
            );
            rightEigenVectors[k] += scaledBase;
          }

          // compute residuum
          V residuum( h.rightApply(rightEigenVectors[k]) );
          V lambdaR(
            rightEigenVectors[k] * ComplexTraits<F>::convert(eigenValues[k])
          );
          residuum -= lambdaR;
          rms += std::real(residuum.dot(residuum)) /
            std::real(rightEigenVectors[k].dot(rightEigenVectors[k]));

          // compute correction using preconditioner
          V correction( p.getCorrection(eigenValues[k], residuum) );

          // orthonormalize and append to basis
          for (unsigned int b(0); b < basis.size(); ++b) {
            V scaledBase( basis[b] * basis[b].dot(correction) );
            correction -= scaledBase;
          }
          correction *= F(1) / std::sqrt(correction.dot(correction));
          basis.push_back(correction);
        }
        ++iterationCount;
        // end basis extension loop
      } while (
        rms >= eigenVectorsCount * tolerance && basis.size() <= maxBasisSize
      );
      // end convergence loop
      if (basis.size() > maxBasisSize) {
        //throw EXCEPTION("Failed to reach convergence");
      }
    }

    const std::vector<complex> &getEigenValues() const {
      return eigenValues;
    }

    const std::vector<V> &getRightEigenVectors() const {
      return rightEigenVectors;
    }

  protected:
    std::vector<complex> eigenValues;
    std::vector<V> rightEigenVectors;
  };
}
#endif

/*
  MySimpleVector v;
  MySimpleMatrix m;
  MySimplePreconditioner<MySimpleMatrix> p
  EigenSystemDavidson<MyFuckingVector> eigenSystem(m, 4, p);
  eigenSystem.getEigenValues()[0];
  eigenSystem.getRightEigenVectors()[0];
*/
