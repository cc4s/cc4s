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
    void testDLogZH0(const unsigned int n = 0, const bool dbeta = D_BETA);

    real beta, deltaMu;
    PTR(CTF::Tensor<>) Dabij, Dai;
  };

  /**
   * \brief Provides a transformation function for the nth derivative of the
   * thermal contraction t = t * 1/(1+exp(-/+eps*beta))
   * either d/(-dbeta) or d/(beta*dmu), -/+ is used for particles/holes,
   * respectively.
   **/
  template <typename F=real>
  class ThermalContraction {
  public:
    ThermalContraction(
      const real beta_,
      const real deltaMu_,
      const bool particle,
      const unsigned int n = 0, const bool dbeta_ = true
    ):
      beta(beta_), deltaMu(deltaMu_), sign(particle ? -1 : +1),
      dbeta(dbeta_), a(n)
    {
      if (n > 0) {
        std::vector<int64_t> nextA(n);
        a[0] = 1;
        for (unsigned int m(3); m <= n+1; ++m) {
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
        if (dbeta) {
          // d/dbeta
          t *= std::pow(sign*eps,a.size()) * y;
        } else {
          // d/dmu
          t *= std::pow(sign,a.size()) * y;
        }
      }
    }
  protected:
    real beta, deltaMu;
    int sign;
    bool dbeta;
    std::vector<int64_t> a;
  };
}

#endif

