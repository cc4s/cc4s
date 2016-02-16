/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COULOMB_VERTEX_DECOMPOSITION_DEFINED
#define COULOMB_VERTEX_DECOMPOSITION_DEFINED

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
  class CoulombVertexDecomposition: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CoulombVertexDecomposition);
    CoulombVertexDecomposition(
      std::vector<Argument> const &argumentList
    );
    virtual ~CoulombVertexDecomposition();
    virtual void run();
      
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
     * \brief The full Coulomb vertex \f$\Gamma^q_{rG}\f$.
     */
    CTF::Tensor<complex> *GammaGqr;
    /**
     * \brief The fit \f${\Pi^\ast}^{qR}\Pi_{rR}\Lambda_{GR}\f$.
     */
    CTF::Tensor<complex> *Gamma0Gqr;
    /**
     * \brief The factor orbitals \f${\Pi^\ast}^{qR} = (\Pi_{qR})^\ast\f$.
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
     * \brief Estimators for the regularization parameter during
     * the alternating least squares fits. They estimate the
     * regularization parameter \f$\lambda\f$ in each iteration from
     * the swamping factor in the previous iteration.
     */
    AlternatingLeastSquaresRegularizationEstimator
      *regularizationEstimatorPiqR, *regularizationEstimatorPirR,
      *regularizationEstimatorLambdaGR;

    static int64_t constexpr DEFAULT_MAX_ITERATIONS = 32;
    static double constexpr DEFAULT_DELTA = 0.0;
    static double constexpr DEFAULT_SWAMPING_THRESHOLD = 1.0;
    static double constexpr DEFAULT_REGULARIZATION_FRICTION = 0.125;
    static bool constexpr DEFAULT_REAL_FACTOR_ORBITALS = false;
    static bool constexpr DEFAULT_NORMALIZED_FACTOR_ORBITALS = false;

  protected:
    /**
     * \brief Performs one iteration in fitting the factor orbitals
     * and the Coulomb factors according to the given algorithm.
     */
    void fit(int64_t iterationsCount);
    /**
     * \brief Normalizes the given factor orbitals, such that
     * \f${\Pi^\ast}^{qR}\Pi_{qR} = \delta_{qq}\f$.
     */
    void normalizePi(CTF::Matrix<complex> &Pi);
    /**
     * \brief Discards the imaginary part of the given factor orbitals.
     */
    void realizePi(CTF::Matrix<complex> &Pi);
  };
}

#endif

