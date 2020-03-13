/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef THERMAL_PERTURBATION_DEFINED 
#define THERMAL_PERTURBATION_DEFINED

#include <algorithms/Algorithm.hpp>

#include <util/SharedPointer.hpp>

namespace cc4s {
  class ThermalPerturbation: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ThermalPerturbation);
    ThermalPerturbation(
      std::vector<Argument> const &argumentList
    );
    virtual ~ThermalPerturbation();
    virtual void run();
    virtual void dryRun();

  protected:
    real getDLogZH0(
      const unsigned int dbeta_n = 0, const unsigned int dmu_m = 0
    );

    real beta, deltaMu;
    PTR(CTF::Tensor<>) Dabij, Dai;
  };

  /**
   * \brief Provides a transformation function for the n/mth derivative of the
   * thermal contraction t = t * 1/(1+exp(-/+eps*beta))
   * with respect to d^n/(-dbeta)^n and d^m/(beta*dmu)^m, respectivel,
   * -/+ is used for particles/holes, respectively.
   **/
  template <typename F=real>
  class ThermalContraction {
  public:
    ThermalContraction(
      const real beta_,
      const real deltaMu_,
      const bool particle,
      const unsigned int dbeta_n_ = 0, const bool dmu_m_ = 0
    ):
      beta(beta_), deltaMu(deltaMu_), sign(particle ? -1 : +1),
      dbeta_n(dbeta_n_), dmu_m(dmu_m_), a(dbeta_n_+dmu_m_)
    {
      if (a.size() > 0) {
        std::vector<int64_t> nextA(a.size());
        a[0] = 1;
        for (unsigned int m(3); m <= a.size()+1; ++m) {
          nextA[0] = 1;
          for (unsigned int k(1); k < m-1; ++k) {
            // build coefficients for nth derivative
            nextA[k] = (k+1)*a[k] - (m-k-1)*a[k-1];
          }
          a = nextA;
        }
      }
    }
    void operator()(const real epsilon, F &t) {
      const real eps(epsilon-deltaMu);
      F x( 1/(1+std::exp(sign*eps*beta)) );
      if (a.size() == 0) {
        // 0th derivative
        t *= x;
      } else {
        // 1st or higher derivative
        F y(0);
        for (size_t k(0); k < a.size(); ++k) {
          // TODO: optimize
          // compose from precomputed coefficients
          y += a[k] * std::pow(x,k+1) * std::pow(1-x,a.size()-k);
        }
        t *= y * std::pow(sign*eps,dbeta_n) * std::pow(sign,dmu_m);
      }
    }
  protected:
    real beta, deltaMu;
    int sign, dbeta_n, dmu_m;
    std::vector<int64_t> a;
  };
}

#endif

