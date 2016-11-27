/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COULOMB_VERTEX_DECOMPOSITION_DEFINED
#define COULOMB_VERTEX_DECOMPOSITION_DEFINED

#include <algorithms/Algorithm.hpp>
#include <tcc/DryTensor.hpp>
#include <math/Complex.hpp>
#include <math/RegularizedAlternatingLeastSquares.hpp>
#include <ctf.hpp>

namespace cc4s {
  /**
   * \brief Decomposes the Coulomb vertex \f$\Gamma_{rG}^q\f$
   * into the factor orbitals \f$\Pi_{qR}\f$
   * and the Coulom factors \f$\Lambda_{GR}\f$.
   * The decomposition is done with a
   * regularized alternating least squares (RALS) algorithm, requiring
   * only a few dozen steps for sufficient convergence.
   * Note that currently the employed fit is \f$\Pi_{qR}\Pi_{rR}\Lambda_{GR}\f$
   * rather than the form with \f$\Pi_{qR}\f$ conjugated.
   */
  class CoulombVertexDecomposition: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CoulombVertexDecomposition);
    CoulombVertexDecomposition(
      std::vector<Argument> const &argumentList
    );
    virtual ~CoulombVertexDecomposition();
    virtual void run();
    virtual void dryRun();
      
    /**
     * \brief The rank \f$N_R\f$ of the tensor rank decomposition
     */
    int64_t rank;

    static int64_t constexpr DEFAULT_RANK = -1;

    static double constexpr DEFAULT_RANK_FACTOR = 2.0;

    /**
     * \brief The Frobenius norm of the difference between
     * \f$\Gamma^q_{rG}\f$ and its decomposition.
     */
    double Delta;

    /**
     * \brief Whether the factor orbitals \f$\Pi_{qR}\f$
     * are required to be real.
     */
    bool realFactorOrbitals;
    /**
     * \brief Whether the factor orbitals \f$\Pi_{qR}\f$
     * are required to be normalized, i.e.
     * \f${\Pi^\ast}^{qR}\Pi_{qR} = \delta{qq}\f$.
     */
    bool normalizedFactorOrbitals;

    /**
     * \brief Whether to write Delta after each part
     * of one RALS iteration, where Delta is the
     * Frobenius norm of
     * \f${\Pi^\ast}^{qR}Pi_{rR}\Lambda_{GR} - \Gamma^q_{rG}\f$.
     */
    bool writeSubIterations;

    /**
     * \brief The full Coulomb vertex \f$\Gamma^q_{rG}\f$.
     */
    CTF::Tensor<complex> *GammaGqr;
    /**
     * \brief The fit \f${\Pi^\ast}^{qR}\Pi_{rR}\Lambda_{GR}\f$.
     */
    CTF::Tensor<complex> *composedGammaGqr;
    /**
     * \brief The conjugated factor orbitals
     * \f${\Pi^\ast}^{qR} = (\Pi_{qR})^\ast\f$.
     */
    CTF::Matrix<complex> *PiqR;
    /**
     * \brief The factor orbitals \f$\Pi_{rR} = ({\Pi^\ast}^{rR})^\ast\f$.
     */
    CTF::Matrix<complex> *PirR;
    /**
     * \brief The full Coulomb factors \f$\Lambda_{GR}\f$.
     */
    CTF::Matrix<complex> *LambdaGR;

    /**
     * \brief Estimator for the regularization parameter during
     * the alternating least squares fits. They estimate the
     * regularization parameter \f$\lambda\f$ in each iteration from
     * the swamping factor in the previous iteration.
     */
    AlternatingLeastSquaresRegularizationEstimator *regularizationEstimator;

    static int64_t constexpr DEFAULT_MAX_ITERATIONS = 32;
    static double constexpr DEFAULT_DELTA = 0.0;
    static double constexpr DEFAULT_SWAMPING_THRESHOLD = 1.0;
    static double constexpr DEFAULT_REGULARIZATION_FRICTION = 0.125;
    static bool constexpr DEFAULT_REAL_FACTOR_ORBITALS = false;
    static bool constexpr DEFAULT_NORMALIZED_FACTOR_ORBITALS = false;
    static bool constexpr DEFAULT_WRITE_SUB_ITERATIONS = false;

  protected:
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
      DryTensor<complex> *GammaGqr,
      DryTensor<complex> *PiqR, DryTensor<complex> *PirR,
      DryTensor<complex> *LambdaGR,
      DryTensor<complex> *composedGammaGqr
    );
    /**
     * \brief Normalizes the given factor orbitals, such that
     * \f${\Pi^\ast}^{qR}\Pi_{qR} = \delta_{qq}\f$.
     */
    void normalizePi(CTF::Matrix<complex> &Pi);
    /**
     * \brief Discards the imaginary part of the given factor orbitals.
     */
    void realizePi(CTF::Matrix<complex> &Pi);

    /**
     * \brief Solves the quadratically occurring factor Pi iteratively
     * similar to the Babylonian algorithm.
     * \note{
     *   Currently \f$\Pi_{qR}Pi_{rR}\Lambda_{GR}\f$ is solved for
     *   instead of \f${\Pi^\ast}^{qR}Pi_{rR}\Lambda_{GR}\f$.
     * }
     */
    void iterateQuadraticFactor(int iterationsCount);
    double getDelta();
  };
}

#endif

