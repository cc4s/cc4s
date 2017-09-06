/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LAPLACE_MP2_ENERGY_DEFINED
#define LAPLACE_MP2_ENERGY_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/Complex.hpp>
#include <tcc/Tensor.hpp>

namespace cc4s {
  /**
   * \brief Caclulates MP2 energy from the Coulomb Integrals \f$V_{ij}^{ab}.
   */
  class LaplaceMp2Energy: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(LaplaceMp2Energy);
    LaplaceMp2Energy(std::vector<Argument> const &argumentList);
    virtual ~LaplaceMp2Energy();
    /**
     * \brief Calculates MP2 energy from Coulomb integrals Vabij
     */
    virtual void run();
    /**
     * \brief Dry run for the MP2 energy from Coulomb integrals Vabij
     */
    virtual void dryRun();

  protected:
    CTF::Tensor<complex> *GpRSn, *GhRSn, *VRS, *wn;
    complex *Gp, *Gh, *V, *w;
    int NR, Nn;

    template <typename F>
    class SampledDistribution {
    public:
      SampledDistribution(
        F *values, int64_t valuesCount_
      ):
        valuesCount(valuesCount_),
        cummulants(new double[valuesCount_])
      {
        double cummulant(0);
        for (int64_t i(0); i < valuesCount_; ++i) {
          cummulants[i] = cummulant;
          cummulant += std::abs(values[i]);
        }
        if (cummulant > 1e-16) {
          for (int64_t i(0); i < valuesCount_; ++i) {
            cummulants[i] /= cummulant;
          }
        }
      }
      ~SampledDistribution() {
        delete[] cummulants;
      }
      void draw(double uniform, int64_t &index, double &weight) {
        int64_t first(0), last(valuesCount-1);
        while (first < last) {
          int64_t mid((first + last + 1) >> 1);
          if (cummulants[mid] <= uniform) {
            first = mid;
          } else {
            last = mid - 1;
          }
        }
        index = first;
        if (first+1 < valuesCount) {
          if ((cummulants[first+1] - cummulants[first]) < 1e-16) {
            throw new EXCEPTION("zero probability");
          }
          weight = 1 / (cummulants[first+1] - cummulants[first]);
        } else {
          weight = 1 / (1 - cummulants[first]);
          if ((1 - cummulants[first]) < 1e-16) {
            throw new EXCEPTION("zero probability");
          }
        }
      }
      int64_t valuesCount;
      double *cummulants;
    };

    void normalizeV(CTF::Tensor<complex> &Lambda, CTF::Tensor<complex> &Pi);
    double calculateNumerically();
    double calculateAnalytically();
    double calculateStochastically();
    double sumNaively();
    double sumMonteCarlo();
    complex getPermutedSamples(int R, int S, int T, int U);
    complex getIntegratedSamples(int R, int S, int T, int U);
    complex getSample(int R, int S, int T, int U, int n);
  };
}

#endif

