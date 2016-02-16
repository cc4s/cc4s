/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef RALS_COULOMB_VERTEX_DECOMPOSITION_DEFINED
#define RALS_COULOMB_VERTEX_DECOMPOSITION_DEFINED

#include <Algorithm.hpp>
#include <math/Complex.hpp>
#include <ctf.hpp>

namespace cc4s {
  /**
   * \brief Decomposes the particle hole Coulomb vertex \f$\Gamma_{iG}^a\f$
   * into the occupied factor orbitals \f$\Pi_{iR}\f$,
   * the virtual factor orbitals \f$\Pi_{aR}\f$, and the Coulom
   * factors \f$\Lambda_{GR}\f$. The decomposition is done with a
   * regularized alternating least squares (RALS) algorithm, requiring
   * only a few dozen steps for sufficient convergence.
   */
  class RalsCoulombVertexDecomposition: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(RalsCoulombVertexDecomposition);
    RalsCoulombVertexDecomposition(
      std::vector<Argument> const &argumentList
    );
    virtual ~RalsCoulombVertexDecomposition();
    virtual void run();
      
    /**
     * \brief The rank of the tensor rank decomposition
     */
    int64_t rank;

    /**
     * \brief The sum of squares of the difference between the current
     * decomposition and the particle hole Coulomb vertex.
     */
    double residuum;
    CTF::Tensor<complex> *GammaGai, *Gamma0Gai;
    CTF::Matrix<complex> *PiiR, *PiaR, *LambdaGR;

    class RegularizationEstimator {
    public:
      RegularizationEstimator(
        double swampingThreshold_, double regularizationFriction_,
        double initialLambda_
      ):
        swampingThreshold(swampingThreshold_),
        regularizationFriction(regularizationFriction_),
        lambda(initialLambda_)
      { }
      double getSwampingThreshold() {
        return swampingThreshold;
      }
      double getLambda() {
        return lambda;
      }
      void update(double const swampingFactor) {
        double s(swampingFactor / swampingThreshold);
        double estimatedLambda(lambda * s*s);
        lambda =
          (1-regularizationFriction)*estimatedLambda +
          regularizationFriction*lambda;
      }
    protected:
      double swampingThreshold, regularizationFriction;
      double lambda;
    };

    RegularizationEstimator
      *regularizationEstimatorPiiR, *regularizationEstimatorPiaR,
      *regularizationEstimatorLambdaGR;

    static int64_t constexpr DEFAULT_MAX_ITERATIONS = 32;
    static double constexpr DEFAULT_DELTA = 0.0;
    static double constexpr DEFAULT_SWAMPING_THRESHOLD = 1.0;
    static double constexpr DEFAULT_REGULARIZATION_FRICTION = 0.125;

  protected:
    void fit(int64_t iterationsCount);
    void calculateGamma0();

    void fitAls(
      char const *indicesGamma,
      CTF::Tensor<complex> &B, char const idxB,
      CTF::Tensor<complex> &C, char const idxC,
      CTF::Tensor<complex> &A, char const idxA
    );
    void fitRals(
      char const *indicesGamma,
      CTF::Tensor<complex> &B, char const idxB,
      CTF::Tensor<complex> &C, char const idxC,
      CTF::Tensor<complex> &A, char const idxA,
      RegularizationEstimator &regularizationEstimatorA
    );

    void applyToGamma(
      char const *indicesGamma,
      CTF::Tensor<complex> &conjB, char const idxB,
      CTF::Tensor<complex> &conjC, char const idxC,
      CTF::Tensor<complex> &A, char const idxA
    );
    void applyToGammaSliced(
      char const *indicesGamma,
      CTF::Tensor<complex> &conjB, char const idxB,
      CTF::Tensor<complex> &conjC, char const idxC,
      CTF::Tensor<complex> &A, char const idxA
    );

    void normalize(CTF::Matrix<complex> &Pi);
    void realize(CTF::Matrix<complex> &Pi);
  };
}

#endif

