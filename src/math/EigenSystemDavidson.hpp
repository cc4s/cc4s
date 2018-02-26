#ifndef EIGEN_SYSTEM_DAVIDSON_DEFINED
#define EIGEN_SYSTEM_DAVIDSON_DEFINED

#include <util/LapackMatrix.hpp>
#include <util/LapackGeneralEigenSystem.hpp>
#include <math/MathFunctions.hpp>

#include <vector>
#include <iomanip>
#include <utility>
#include <algorithm>

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
     * Also, if the dual version of the Davidson algorithm is to be used
     * h should offer the leftApply method
     * V leftApply(V &v);
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
     * \param[in] dualVersion If the dual version of the algorithm is to be
     * used. The dual version of the algorithm calculates both the right and
     * left eigenvectors of the h object. If the dual version is not to be used
     * then only the right eigenvectors are to be calculated.
     * \param[in] refreshIterations This vector of integers represents
     * the iterations where the refreshment should be done. Refreshment meaning
     * that the trial Davidson basis should be thrown away and only keep
     * the current approximations to the eigenvectors.
     **/
    template <typename H, typename P>
    EigenSystemDavidson(
      H &h,
      const int eigenVectorsCount,
      P &p,
      const double tolerance = 1E-14,
      const unsigned int maxBasisSize = 1000,
      const unsigned int maxIterations = 1000,
      const unsigned int minIterations = 1,
      const bool dualVersion = false,
      std::vector<int> refreshIterations = std::vector<int>{{}}
    ):
      eigenValues(eigenVectorsCount)
    {
      if (dualVersion)
        eigenSystemDualVersion(
          h, eigenVectorsCount, p, tolerance, maxBasisSize,
          maxIterations, minIterations
        );
      else
        eigenSystemMonoVersion(
          h, eigenVectorsCount, p, tolerance, maxBasisSize,
          maxIterations, minIterations, refreshIterations

        );
    }

    template <typename H, typename P>
    void eigenSystemDualVersion(
      H &h,
      const int eigenVectorsCount,
      P &p,
      const double tolerance,
      const unsigned int maxBasisSize,
      const unsigned int maxIterations,
      const unsigned int minIterations
    ) {
      // get inital estimates for rEV = initial B matrix
      rightEigenVectors = p.getInitialBasis(eigenVectorsCount);
      leftEigenVectors = p.getInitialBasis(eigenVectorsCount);
      std::vector<V> rightBasis( rightEigenVectors );
      std::vector<V> leftBasis( leftEigenVectors );

      // begin convergence loop
      double rms;
      unsigned int iterationCount(0);
      do {
        LOG(1,"Davidson") << "iteration=" << (iterationCount+1) << std::endl;

        // check if leftBasis and rightBasis have the same size
        if (rightBasis.size() != leftBasis.size()) {
          LOG(1, "Davidson") << "Right and left bases have different sizes"
                             << std::endl;
          throw EXCEPTION("");
        }

        // compute reduced H by projection onto subspace spanned by rightBasis
        // and leftBasis
        LapackMatrix<complex> reducedH(rightBasis.size(), rightBasis.size());
        for (unsigned int j(0); j < rightBasis.size(); ++j) {
          V HBj( h.rightApply(rightBasis[j]) );
          for (unsigned int i(0); i < rightBasis.size(); ++i) {
            reducedH(i,j) = leftBasis[i].dot(HBj);
          }
        }

        // compute K lowest reduced eigenvalues and vectors of reduced H
        // also compute left eigenvectors
        LapackGeneralEigenSystem<complex> reducedEigenSystem(
          reducedH, true, true
        );

        // begin rightBasis extension loop for each k
        rms = 0.0;
        for (unsigned int k(0); k < eigenValues.size(); ++k) {
          // get estimated eigenvalue
          eigenValues[k] = reducedEigenSystem.getEigenValues()[k];

          // compute estimated eigenvector by expansion in rightBasis
          rightEigenVectors[k] *= F(0);
          leftEigenVectors[k] *= F(0);
          for (int b(0); b < reducedH.getColumns(); ++b) {
            // Get the rightEigenVectors in the coordinates of the
            // original basis (rightBasis)
            V rightScaledBase(
              rightBasis[b] * Cast<F>::from(
                reducedEigenSystem.getRightEigenVectors()(b,k)
              )
            );
            rightEigenVectors[k] += rightScaledBase;

            // Do the same for the leftEigenVectors
            V leftScaledBase(
              leftBasis[b] * Cast<F>::from(
                reducedEigenSystem.getLeftEigenVectors()(b,k)
              )
            );
            leftEigenVectors[k] += leftScaledBase;

          }

          double rightNorm(
            rightEigenVectors[k].dot(rightEigenVectors[k])
          );
          LOG(1,"Davidson") << "Right norm [" << k << "] "
                            << rightNorm << std::endl;
          LOG(1,"Davidson") << "EV         [" << k << "] "
                            << eigenValues[k] << std::endl;

          // compute residuum
          V rightResiduum( h.rightApply(rightEigenVectors[k]) );
          V rightLambdaR(
            rightEigenVectors[k] * Cast<F>::from(eigenValues[k])
          );
          rightResiduum -= rightLambdaR;

          std::cout << "Applying left" << std::endl;
          V leftResiduum( h.leftApply(leftEigenVectors[k]) );
          V leftLambdaR(
            leftEigenVectors[k] * Cast<F>::from(eigenValues[k])
          );
          leftResiduum -= leftLambdaR;

          rms += std::real(rightResiduum.dot(rightResiduum)) /
            std::real(rightEigenVectors[k].dot(rightEigenVectors[k]));

          // compute rightCorrection using preconditioner
          V rightCorrection( p.getCorrection(eigenValues[k], rightResiduum) );

          // orthonormalize and append to rightBasis
          for (unsigned int b(0); b < rightBasis.size(); ++b) {
            V rightScaledBase( rightBasis[b] * rightBasis[b].dot(rightCorrection) );
            rightCorrection -= rightScaledBase;
          }
          F correction_norm(rightCorrection.dot(rightCorrection));
          if (std::abs(correction_norm) < tolerance) continue;
          rightCorrection *= F(1) / correction_norm;
          rightBasis.push_back(rightCorrection);
        }
        ++iterationCount;
        // end rightBasis extension loop
      } while (
        rms >= eigenVectorsCount * tolerance &&
        rightBasis.size() <= maxBasisSize
      );
      // end convergence loop
      if (rightBasis.size() > maxBasisSize) {
        //throw EXCEPTION("Failed to reach convergence");
      }
    }

    template <typename H, typename P>
    void eigenSystemMonoVersion(
      H &h,
      const int eigenVectorsCount,
      P &p,
      const double tolerance,
      const unsigned int maxBasisSize,
      const unsigned int maxIterations,
      const unsigned int minIterations,
      std::vector<int> refreshIterations
    ) {
      LOG(1,"Davidson").flags(
        std::ios::right | std::ios::scientific | std::ios::showpos
      );
      // get inital estimates for rEV = initial B matrix
      rightEigenVectors = p.getInitialBasis(eigenVectorsCount);

      std::vector<V> rightBasis( rightEigenVectors );

      // begin convergence loop
      double rms;
      unsigned int iterationCount(0);
      do {
        LOG(1,"Davidson") << "iteration=" << (iterationCount+1) << std::endl;

        // Check if a refreshment should be done
        if (
          std::find(
            refreshIterations.begin(),
            refreshIterations.end(),
            iterationCount + 1
          ) != refreshIterations.end()
        ) {
          LOG(1,"Davidson") << "Refreshing BASIS" << std::endl;
          rightBasis.resize(eigenVectorsCount);
          eigenValues.resize(eigenVectorsCount);
          for (unsigned int i(0) ; i < rightBasis.size() ; i++) {
            rightBasis[i] = rightEigenVectors[i];
          }
        }

        // compute reduced H by projection onto subspace spanned by rightBasis
        LapackMatrix<complex> reducedH(rightBasis.size(), rightBasis.size());

        for (unsigned int j(0); j < rightBasis.size(); ++j) {
          V HBj( h.rightApply(rightBasis[j]) );
          for (unsigned int i(0); i < rightBasis.size(); ++i) {
            //V HBi( h.rightApply(rightBasis[i]) );
            //reducedH(i,j) = HBi.dot(HBj);
            reducedH(i,j) = rightBasis[i].dot(HBj);
          }
        }

        // compute K lowest reduced eigenvalues and vectors of reduced H
        LapackMatrix<complex> reducedEigenVectors(
          rightBasis.size(), rightBasis.size()
        );
        LapackGeneralEigenSystem<complex> reducedEigenSystem(reducedH);

        // begin rightBasis extension loop for each k
        rms = 0.0;
        for (unsigned int k(0); k < eigenValues.size(); ++k) {
          // get estimated eigenvalue
          eigenValues[k] = reducedEigenSystem.getEigenValues()[k];

          // compute estimated eigenvector by expansion in rightBasis
          rightEigenVectors[k] *= F(0);
          for (int b(0); b < reducedH.getColumns(); ++b) {
            rightEigenVectors[k] +=
              rightBasis[b] * Cast<F>::from(
                reducedEigenSystem.getRightEigenVectors()(b,k)
              );
          }
          double rightNorm(
            rightEigenVectors[k].dot(rightEigenVectors[k])
          );

          LOG(1,"Davidson") << "Right norm [" << k << "] "
                            << rightNorm << std::endl;
          LOG(1,"Davidson") << "EV         [" << k << "] "
                            << std::setprecision(15) << std::setw(23)
                            << eigenValues[k] << std::endl;



          // compute residuum
          V residuum( h.rightApply(rightEigenVectors[k]) );
          residuum -=
            rightEigenVectors[k] * Cast<F>::from(
              //std::sqrt(eigenValues[k])
              eigenValues[k]
            );
          rms += std::real(residuum.dot(residuum)) /
            std::real(rightEigenVectors[k].dot(rightEigenVectors[k]));

          // compute correction using preconditioner
          V correction( p.getCorrection(eigenValues[k], residuum) );

          // orthonormalize and append to rightBasis
          for (unsigned int b(0); b < rightBasis.size(); ++b) {
            correction -= rightBasis[b] * rightBasis[b].dot(correction);
          }
          F correction_norm( std::sqrt(correction.dot(correction)) );
          if (std::abs(correction_norm) < 1E-6) continue;
          correction *= 1 / correction_norm;
          rightBasis.push_back(correction);
          LOG(1,"Davidson") << "Basis size " << rightBasis.size() << std::endl;
        }
        ++iterationCount;
        // end rightBasis extension loop
      } while (
        iterationCount+1 <= minIterations    || (
          rms >= eigenVectorsCount * tolerance &&
          rightBasis.size() <= maxBasisSize    &&
          iterationCount+1 <= maxIterations
        )
      );
      // end convergence loop
      if (rightBasis.size() > maxBasisSize) {
        //throw EXCEPTION("Failed to reach convergence");
      }
    }

    const std::vector<complex> &getEigenValues() const {
      return eigenValues;
    }

    const std::vector<V> &getRightEigenVectors() const {
      return rightEigenVectors;
    }

    const std::vector<V> &getLeftEigenVectors() const {
      return leftEigenVectors;
    }

  protected:
    std::vector<complex> eigenValues;
    std::vector<V> rightEigenVectors;
    std::vector<V> leftEigenVectors;
  };
}
#endif

/*
  MySimpleVector v;
  MySimpleMatrix m;
  MySimplePreconditioner<MySimpleMatrix> p
  EigenSystemDavidson<MyVectorType> eigenSystem(m, 4, p);
  eigenSystem.getEigenValues()[0];
  eigenSystem.getRightEigenVectors()[0];
*/
