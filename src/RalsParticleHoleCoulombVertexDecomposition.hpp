/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef RALS_PARTICLE_HOLE_COULOMB_VERTEX_DECOMPOSITION_DEFINED
#define RALS_PARTICLE_HOLE_COULOMB_VERTEX_DECOMPOSITION_DEFINED

#include <Algorithm.hpp>
#include <util/Complex.hpp>
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
    CTF::Tensor<complex> *gammaGai, *gamma0Gai;
    CTF::Matrix<complex> *piiR, *piaR, *lambdaGR;

  protected:
    void fit(double lambda);

    void fitAls(
      char const *indicesGamma,
      CTF::Tensor<complex> &b, char const idxB,
      CTF::Tensor<complex> &c, char const idxC,
      CTF::Tensor<complex> &a, char const idxA
    );
    double fitRals(
      char const *indicesGamma,
      CTF::Tensor<complex> &b, char const idxB,
      CTF::Tensor<complex> &c, char const idxC,
      CTF::Tensor<complex> &a, char const idxA,
      double lambda
    );

    void normalizePi(CTF::Matrix<complex> &pi);
    void realizePi(CTF::Matrix<complex> &pi);
  };
}

#endif

