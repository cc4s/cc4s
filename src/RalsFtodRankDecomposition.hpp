/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef RALS_FTOD_RANK_DECOMPOSITION_DEFINED
#define RALS_FTOD_RANK_DECOMPOSITION_DEFINED

#include <Algorithm.hpp>
#include <util/Complex.hpp>
#include <ctf.hpp>


namespace cc4s {
  /**
   * \brief This algorithm provides a tensor rank decomposition of the
   * Fourier tranformed overlap densities.
   */
  class RalsFtodRankDecomposition: public Algorithm {
  public:
    RalsFtodRankDecomposition(std::vector<Argument const *> const &argumentList);
    virtual ~RalsFtodRankDecomposition();
    /**
     * \deprecated
     */
    virtual std::vector<std::string> getDefaultArgumentOrder() {
      std::vector<std::string> argumentOrder;
      return argumentOrder;
    }
    virtual void run();
      
    /**
     * \brief the rank of the tensor rank decomposition
     */
    int64_t rank;
    double R;
    CTF::Tensor<complex> *chi, *chi0;
    CTF::Matrix<complex> *x, *gamma;

    static void test(CTF::World *world);
  protected:
    void fit(double lambda);

    void fitAls(
      char const *indicesChi,
      CTF::Tensor<complex> &b, char const idxB,
      CTF::Tensor<complex> &c, char const idxC,
      CTF::Tensor<complex> &a, char const idxA
    );
    double fitRals(
      char const *indicesChi,
      CTF::Tensor<complex> &b, char const idxB,
      CTF::Tensor<complex> &c, char const idxC,
      CTF::Tensor<complex> &a, char const idxA,
      double lambda
    );
  };
}

#endif

