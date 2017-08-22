#ifndef EIGEN_SYSTEM_DAVIDSON_DEFINED
#define EIGEN_SYSTEM_DAVIDSON_DEFINED

#include <util/LapackMatrix.hpp>
#include <util/LapackGeneralEigenSystem.hpp>

#include <vector>

namespace cc4s {
  template <typename V>
  class EigenSystemDavidson {
  public:
    typedef typename V::FieldType F;

    /**
     * \brief ...
     * \param[in] h object representing the matrix whose eigen system is sought
     * offering the following method:
     * V rightApply(const V &v) const;
     * returning the action of h on v standing right of h.
     * \param[in] eigenVectorsCount the number of eigenvalues and vectors to be
     * computed.
     * \param[in] p object representing the preconditioner
     * offering the following method:
     * std::vector<V> getInitialBasis(int eigenVectorsCount) const;
     * returning an initial guess the basis consisting of eigenVectorsCount
     * vectors.
     * V getCorrection(const complex eigenValue, const V &residuum) const;
     * computing the estimated correction for the k-th eigenvector given its
     * eigenvalue and the residuum of the current k-th estimated eigenvector.
     * \param[in] tolerance the targeted relative tolerance to be met
     * by all residua.
     * \param[in] maxBasisSize the maximum allowed number of vectors
     * representing eigenvectors.
     **/
    template <typename H, typename P>
    EigenSystemDavidson(
      const H &h,
      const int eigenVectorsCount,
      const P &p,
      const double tolerance = 1E-14,
      const unsigned int maxBasisSize = 1000
    ) {
      eigenValues.resize(eigenVectorsCount);
      rightEigenVectors.resize(eigenVectorsCount);

      // get inital estimates for rEV = initial B matrix
      std::vector<V> basis;
      basis = p.getInitialBasis(eigenVectorsCount);

      // begin convergence loop
      double rms;
      do {
        // compute reduced H by projection onto subspace spanned by basis
        LapackMatrix<F> reducedH(basis.size(), basis.size());
        for (unsigned int j(0); j < basis.size(); ++j) {
          V HBj( h.rightApply(basis[j]) );
          for (unsigned int i(0); i < basis.size(); ++i) {
            reducedH(i,j) = basis[i].dot(HBj);
          }
        }

        // compute K lowest reduced eigenvalues and vectors of reduced H
        LapackMatrix<F> reducedEigenVectors(basis.size(), basis.size());
        LapackGeneralEigenSystem<F> reducedEigenSystem(reducedH);
        
        // begin basis extension loop for each k
        rms = 0.0;
        for (unsigned int k(0); k < eigenValues.size(); ++k) {
          // get estimated eigenvalue
          F estimatedEigenValue( reducedEigenSystem.getEigenValues()[k] );

          // compute estimated eigenvector by expansion in basis
          V estimatedEigenVector(rightEigenVectors[0]);
          estimatedEigenVector *= F(0);
          for (int b(0); b < reducedH.getColumns(); ++b) {
            estimatedEigenVector +=
              basis[b] * reducedEigenSystem.getRightEigenVectors()(b,k);
          }

          // compute residuum
          V residuum( h.rightApply(estimatedEigenVector) );
          residuum -= estimatedEigenVector * estimatedEigenValue;
          rms += std::real(residuum.dot(residuum)) /
            std::real(estimatedEigenVector.dot(estimatedEigenVector));

          // compute correction using preconditioner
          V correction( p.getCorrection(estimatedEigenValue, residuum) );

          // orthonormalize and append to basis
          for (unsigned int b(0); b < basis.size(); ++b) {
            correction -= basis[b] * basis[b].dot(correction);
          }
          correction *= F(1) / std::sqrt(correction.dot(correction));
          basis.push_back(correction);
        }
        // end basis extension loop
      } while (
        rms >= eigenVectorsCount * tolerance && basis.size() <= maxBasisSize
      );
      // end convergence loop
      if (basis.size() > maxBasisSize) {
        throw EXCEPTION("Failed to reach convergence");
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
