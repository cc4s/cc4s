/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PARTICLE_HOLE_COULOMB_VERTEX_DECOMPOSITION_DEFINED
#define PARTICLE_HOLE_COULOMB_VERTEX_DECOMPOSITION_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/Complex.hpp>
#include <math/RegularizedAlternatingLeastSquares.hpp>
#include <tcc/DryTensor.hpp>
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
  class ParticleHoleCoulombVertexDecomposition: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ParticleHoleCoulombVertexDecomposition);
    ParticleHoleCoulombVertexDecomposition(
      std::vector<Argument> const &argumentList
    );
    virtual ~ParticleHoleCoulombVertexDecomposition();
    virtual void run();
    virtual void dryRun();
      
    /**
     * \brief The rank \f$N_R\f$ of the tensor rank decomposition
     */
    int64_t rank;

    /**
     * \brief The Frobenius norm of the difference between
     * \f$\Gamma^a_{iG}\f$ and its decomposition.
     */
    double Delta;

    /**
     * \brief The number of steps between two evaluations of the
     * error in the MP2 energy of the decomposition.
     * A value below 1 results in no evaluations.
     */
    int epsilonStep;

    /**
     * \brief Whether the factor orbitals \f$\Pi_{aR},\Pi_{iR}\f$
     * are required to be real.
     */
    bool realFactorOrbitals;
    /**
     * \brief Whether the factor orbitals \f$\Pi_{aR},\Pi_{iR}\f$
     * are required to be normalized, i.e.
     * \f${\Pi^\ast}^{qR}\Pi_{qR} = \delta{qq}\f$.
     */
    bool normalizedFactorOrbitals;

    /**
     * \brief The Coulomb vertex \f$\Gamma^a_{iG}\f$ restricted to
     * a particle and a hole index.
     */
    CTF::Tensor<complex> *GammaGai;
    /**
     * \brief The fit \f${\Pi^\ast}^{aR}\Pi_{iR}\Lambda_{GR}\f$.
     */
    CTF::Tensor<complex> *Gamma0Gai;
    /**
     * \brief The factor orbitals for occupied states \f$\Pi_{iR}\f$.
     */
    CTF::Tensor<complex> *PiiR;
    /**
     * \brief The conjugated factor orbitals for virtual states
     * \f${\Pi^\ast}^{aR}\f$.
     */
    CTF::Tensor<complex> *PiaR;
    /**
     * \brief The Coulomb factors \f$\Lambda_{GR}\f$
     * in the particle/hole decomposition.
     */
    CTF::Tensor<complex> *LambdaGR;

    /**
     * \brief Estimators for the regularization parameter during
     * the alternating least squares fits. They estimate the
     * regularization parameter \f$\lambda\f$ in each iteration from
     * the swamping factor in the previous iteration.
     */
    AlternatingLeastSquaresRegularizationEstimator
      *regularizationEstimatorPiiR, *regularizationEstimatorPiaR,
      *regularizationEstimatorLambdaGR;

    static int64_t constexpr DEFAULT_MAX_ITERATIONS = 32;
    static double constexpr DEFAULT_DELTA = 0.0;
    static int constexpr DEFAULT_EPSILON_STEP = 0;
    static double constexpr DEFAULT_SWAMPING_THRESHOLD = 1.0;
    static double constexpr DEFAULT_REGULARIZATION_FRICTION = 0.125;
    static bool constexpr DEFAULT_REAL_FACTOR_ORBITALS = false;
    static bool constexpr DEFAULT_NORMALIZED_FACTOR_ORBITALS = false;

  protected:
    double mp2Energy;

    /**
     * \brief Performs one iteration in fitting the factor orbitals
     * and the Coulomb factors according to the given algorithm.
     */
    void fit(int64_t iterationsCount);
    /**
     * \brief Performs a dry run of one iteration in fitting the factor
     * orbitals and the Coulomb factors according to the given algorithm.
     */
    void dryFit(
      DryTensor<complex> *GammaGai,
      DryTensor<complex> *PiaR, DryTensor<complex> *PiiR,
      DryTensor<complex> *LambdaGR,
      DryTensor<complex> *Gamma0Gai
    );
    /**
     * \brief Normalizes the given factor orbitals, such that
     * \f${\Pi^\ast}^{qR}\Pi_{qR} = \delta_{qq}\f$.
     */
    void normalizePi(CTF::Tensor<complex> &Pi);
    /**
     * \brief Discards the imaginary part of the given factor orbitals.
     */
    void realizePi(CTF::Tensor<complex> &Pi);

    /**
     * \brief Evaluates and prints the error of the MP2 energy between
     * the current decomposition and the full MP2 energy.
     */
    void evaluateMp2Error();
    /**
     * \brief Evaluates the MP2 energy for the given Coulomb vertex.
     */
    double evaluateMp2(CTF::Tensor<complex> &Gamma);
    /**
     * \brief Performs a dry run of evaluating the MP2 energy
     * for the given Coulomb vertex.
     */
    void dryEvaluateMp2(DryTensor<complex> &Gamma);
  };
}

#endif

