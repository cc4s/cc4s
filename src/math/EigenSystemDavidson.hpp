#ifndef EIGEN_SYSTEM_DAVIDSON_DEFINED
#define EIGEN_SYSTEM_DAVIDSON_DEFINED

#include <util/LapackMatrix.hpp>
#include <util/LapackGeneralEigenSystem.hpp>

#include <vector>

namespace cc4s {
  template <typename V>
  void EigenSystemDavidson {
  public:
    typedef typename V::FieldType F;

    /**
     * \brief ...
     * \param[in] h object representing the matrix whose eigen system is sought
     * offering the following method:
     * V rightApply(const V &v) returning the action of h on v standing right of h.
     * \param[in] eigenVectorsCount the number of eigenvalues and vectors to be
     * computed.
     * \param[in] p object representing the preconditioner
     * offering the following method:
     * std::vector<V> getInitialBasis(int eigenVectorsCount) returning an
     * initial guess the basis consisting of eigenVectorsCount vectors.
     * V correction(const complex eigenValue, const V &residue);
     * computing the estimated correction for the k-th eigenvector given its
     * eigenvalue and the residue of the current k-th estimated eigenvector.
     * \param[in] tolerance the targeted relative tolerance to be met
     * by all residues.
     * \param[in] maxBasisSize the maximum allowed number of vectors
     * representing eigenvectors.
     **/
    template <typename H, typename P>
    EigenSystemDavidson(
      const H &h,
      const int eigenVectorsCount,
      const P &p,
      const double tolerance = 1E-100,
      const int maxBasisSize = 1000
    ) {
      eigenValues.resize(eigenVectorsCount);
      rightEigenVectors.resize(eigenVectorsCount);

      // get inital estimates for rEV = initial B matrix
      std::vector<V> basis;
      basis = p.getInitalBasis(eigenVectorsCount);

      // begin convergence loop
      while (basis.size() <= maxBasisSize) {
        // compute reduced H
        LapackMatrix<F> reducedH(basis.size(), basis.size());
        for (int i(0); i < basis.size(); ++i) {
          for (int j(0); j < basis.size(); ++j) {
            reducedH(i,j) = basis[i].dot( h.rightApply(basis[j]) );
          }
        }
        // compute K lowest reduced eigenvalues and vectors of reduced H
        LapackMatrix<F> reducedEigenVectors(basis.size(), basis.size());
        LapackGeneralEigenSystem<F> reducedEigenSystem(
          reducedH, &reducedEigenVectors
        );
        std::vector<complex> reducedEigenValues(
          reducedEigenSystem.solve()
        );
        
        // begin basis extension loop for each k
        for (int k(0); k < eigenValues.size(); ++k) {
          // compute estimated eigenvector
          V estimatedEigenVector(rightEigenVector[0]);
          estimatedEigenVector *= F(0);
          for (int b(0); b < reducedH.columns; ++b) {
            estimatedEigenVector += basis[b] * reducedEigenVectors(b,k);
          }

          // compute residuum
          V residue(estimatedEigenVector * (-reducedEigenValues[k]));
          residue += h.rightApply(estimatedEigenVector);

          // compute correction using preconditioner
          V correction(p.correction(reducedEigenValues[k], residue));

          // orthonormalize and append to basis
          for (int b(0); b < basis.size(); ++b) {
            correction += basis[b] * ( -basis[b].dot(correction) );
          }
          correction *= F(1) / cc4s::sqrt(correction.dot(correction));
          basis.push_back(correction);
        }
        // end basis extension loop
      }
      // end convergence loop
    }


  protected:
    std::vector<complex> eigenValues;
    std::vector<V> rightEigenVectors;
  };
}
#endif

