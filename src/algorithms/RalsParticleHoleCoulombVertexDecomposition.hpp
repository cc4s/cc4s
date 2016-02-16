/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef RALS_PARTICLE_HOLE_COULOMB_VERTEX_DECOMPOSITION_DEFINED
#define RALS_PARTICLE_HOLE_COULOMB_VERTEX_DECOMPOSITION_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/Complex.hpp>
#include <math/RegularizedAlternatingLeastSquares.hpp>
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
  class RalsParticleHoleCoulombVertexDecomposition: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(RalsParticleHoleCoulombVertexDecomposition);
    RalsParticleHoleCoulombVertexDecomposition(
      std::vector<Argument> const &argumentList
    );
    virtual ~RalsParticleHoleCoulombVertexDecomposition();
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

    AlternatingLeastSquaresRegularizationEstimator
      *regularizationEstimatorPiiR, *regularizationEstimatorPiaR,
      *regularizationEstimatorLambdaGR;

    static int64_t constexpr DEFAULT_MAX_ITERATIONS = 32;
    static double constexpr DEFAULT_DELTA = 0.0;
    static double constexpr DEFAULT_SWAMPING_THRESHOLD = 1.0;
    static double constexpr DEFAULT_REGULARIZATION_FRICTION = 0.125;

  protected:
    void fit(int64_t iterationsCount);

    void fitAls(
      char const *indicesGamma,
      CTF::Tensor<complex> &B, char const idxB,
      CTF::Tensor<complex> &C, char const idxC,
      CTF::Tensor<complex> &A, char const idxA
    );

    void normalizePi(CTF::Matrix<complex> &Pi);
    void realizePi(CTF::Matrix<complex> &Pi);
  };
}

#endif

