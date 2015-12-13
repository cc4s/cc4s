/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CROSS_ENTROPY_FTOD_RANK_DECOMPOSITION_DEFINED
#define CROSS_ENTROPY_FTOD_RANK_DECOMPOSITION_DEFINED

#include <Algorithm.hpp>
#include <random>
#include <ctf.hpp>

namespace cc4s {
  class RankDecomposition {
  public:
    RankDecomposition(int rank, int np, int nG, CTF::World *world);
    CTF::Matrix<> X, GamR, GamI;
    double residuum;
  };

  /**
   * \brief This algorithm provides a tensor rank decomposition of the
   * Fourier tranformed overlap densities.
   */
  class CrossEntropyFtodRankDecomposition: public Algorithm {
  public:
    CrossEntropyFtodRankDecomposition(
      std::vector<Argument> const &argumentList
    );
    virtual ~CrossEntropyFtodRankDecomposition();
    virtual void run();
      
    /**
     * \brief the rank of the tensor rank decomposition
     */
    int64_t rank;
    double R;
    CTF::Tensor<> *chiR, *chiI, *chi0R, *chi0I, *RR, *RI, *XX;
    CTF::Matrix<> *X, *gamR, *gamI;
    std::mt19937 random;

  protected:
    RankDecomposition *estimator, *mu, *sigma;
    int samplesCount, estimatorsCount;
    RankDecomposition **estimators;
    int worstEstimator;

    void calculateChi0(RankDecomposition &d);
    void calculateResiduum();
    double findWorstEstimator();
    void setRandom(
      RankDecomposition &d, RankDecomposition &mu, RankDecomposition &sigma
    );
    void setRandom(CTF::Tensor<> &t, CTF::Tensor<> &mu, CTF::Tensor<> &sigma);
    void calculateMu();
    double calculateSigma();
  };
}

#endif

