/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef FTOD_RANK_DECOMPOSITION_DEFINED
#define FTOD_RANK_DECOMPOSITION_DEFINED

#include <Algorithm.hpp>
#include <ctf.hpp>

namespace cc4s {
  class FtodOptimizationMask {
  public:
    virtual double weight(int q, int r, int G) const = 0;
  };

  /**
   * \brief This algorithm provides a tensor rank decomposition of the
   * Fourier tranformed overlap densities.
   */
  class FtodRankDecomposition: public Algorithm {
  public:
    FtodRankDecomposition(std::vector<Argument const *> const &argumentList);
    virtual ~FtodRankDecomposition();
    virtual void run();
      
    /**
     * \brief the rank of the tensor rank decomposition
     */
    int64_t rank;
    double R;
    CTF::Tensor<> *chiR, *chiI, *chi0R, *chi0I, *RR, *RI;
    CTF::Matrix<> *X, *gamR, *gamI;

  protected:  
    /**
     * \brief gradient of \f$R\f$ with respect to \f$ X\f$.
     */
    CTF::Matrix<> *dX;
    /**
     * \brief gradient of \f$R\f$ with respect to the real part of
     * \f$ \Gamma\f$.
     */
    CTF::Matrix<> *dGamR;
    /**
     * \brief gradient of \f$R\f$ with respect to the imaginary part of
     * \f$ \Gamma\f$.
     */
    CTF::Matrix<> *dGamI;
    /**
     * \brief search direction for \f$X\f$.
     */
    CTF::Matrix<> *sX;
    /**
     * \brief real part of the search direction for \f$\Gamma\f$.
     */
    CTF::Matrix<> *sGamR;
    /**
     * \brief imaginary part of the search direction for \f$\Gamma\f$.
     */
    CTF::Matrix<> *sGamI;

    CTF::Tensor<> *mask;

    /**
     * \brief coefficients of
     * \f${\rm const.} + \alpha^1a_1 + \alpha^2a_2+ \alpha^3a_3+ \alpha^4a_4\f$
     * for the line search of \f$X\f$.
     * \deprecated
     */
    double a0, a1, a2, a3, a4;
    /**
     * \brief coefficients for the line search
     */
    double a[7];

    void setOptimizationMask(FtodOptimizationMask const &mask);
    void normalizeX();
    void initializeRandom(CTF::Tensor<> &t, int64_t seed);
    void initializeX();
    void initializeGam();

    void calculateChi0();
    void calculateResiduum();

    void calculateGradient();

    void lineSearchPart(
      CTF::Tensor<> &R, CTF::Tensor<> &G, CTF::Tensor<> &dG
    );

    double lineSearch();

    void optimize(double const epsilon = 1e-10);

    /**
     * \deprecated
     */
    void lineSearchXPart(
      CTF::Tensor<> &chi0, CTF::Tensor<> &chi, CTF::Tensor<> &gam
    );
    /**
     * \deprecated
     */
    void lineSearchGamPart(
      CTF::Tensor<> &chi0, CTF::Tensor<> &chi, CTF::Tensor<> &dGam
    );

    /**
     * \brief calculates
     * \f${\rm argmin}_\alpha d(X_q(R)+\alpha X_q(R),\Gamma_G(R))\f$
     * where \f$ d(X_q(R),\Gamma_G(R)) \f$ is the (square of) the Euclidian
     * distance betwenn the tensor rank approximation and \f$\chi_r^q(G)\f$.
     * \deprecated
     */
    double lineSearchX();
    /**
     * \brief calculates
     * \f${\rm argmin}_\alpha d(X_q(R),\Gamma_G(R)+\alpha \Gamma_G(R))\f$
     * where \f$ d(X_q(R),\Gamma_G(R)) \f$ is the (square of) the Euclidian
     * distance betwenn the tensor rank approximation and \f$\chi_r^q(G)\f$.
     * \deprecated: rather use the full lineSearch.
     */
    double lineSearchGam();

    void optimizeX(double const epsilon = 1e-10);
    void optimizeGam(double const epsilon = 1e-10);
    // TODO: put in separate test class
    void testGradient();
  };
}

#endif

