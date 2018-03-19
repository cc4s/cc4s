/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef THERMAL_MP2_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED 
#define THERMAL_MP2_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>

#include <util/SharedPointer.hpp>

namespace cc4s {
  class ThermalMp2EnergyFromCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ThermalMp2EnergyFromCoulombIntegrals);
    ThermalMp2EnergyFromCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~ThermalMp2EnergyFromCoulombIntegrals();
    /**
     * \brief Calculates the thermal MP2 energy from Coulomb integrals Vabij
     */
    virtual void run();
    /**
     * \brief Dry run for the thermal MP2 energy from Coulomb integrals Vabij
     */
    virtual void dryRun();

  protected:
    real beta;
    PTR(CTF::Tensor<>) Dabij;

    void computeFreeEnergy();
    void computeEnergyMoments();

    /**
     * \brief Computes the nth derivative of the logarithm of the grand
     * canonical partition function Z(beta) w.r.t. (-beta).
     **/
    real getDerivativeLogZ(const int n = 0);

    real getDerivativeLogZMp2(const int n = 0);
    real getDerivativeLogZHf(const int n = 0);
    real getDerivativeLogZH0(const int n = 0);

    void testDerivativeLogZMp2(const int n = 0);
    void testDerivativeLogZHf(const int n = 0);
    void testDerivativeLogZH0(const int n = 0);

    void addLogZMp2Amplitudes(
      CTF::Tensor<> &Tabij, const std::vector<int> &degrees
    );
    void addLogZHfAmplitudes(
      CTF::Tensor<> &Tij, const std::vector<int> &degrees
    );

    void writeContribution(
      const std::string &contribution, const int n, const real derivativeLogZ
    );
  };

  /**
   * \brief Provides a transformation function for the nth derivative of the
   * thermal second order propagation
   * t = t * d^n/d(-beta)^n int_t1^beta dt2 int_0^beta dt1 exp(-Delta*(t2-t1)).
   **/
  template <typename F=real>
  class ThermalMp2Propagation {
  public:
    ThermalMp2Propagation(
      const real beta_, const int n_ = 0
    ): beta(beta_), n(n_) {
    }
    void operator()(const real Delta, F &t) {
      const real e( std::exp(-beta*Delta) );
      switch (n) {
      case 0:
        if (std::abs(Delta) > 1e-8) {
          t *= (e - 1.0 + Delta*beta) / (Delta*Delta);
        } else {
          t *= beta/2*(beta - beta/3*beta*Delta);
        }
        break;
      case 1:
        if (std::abs(Delta) > 1e-8) {
          t *= (e - 1.0) / Delta;
        } else {
          t *= beta*(-1.0 + beta*Delta/2);
        }
        break;
      default:
        t *= std::pow(Delta,n-2) * e;
      }
    }
  protected:
    real beta;
    int n;
  };

  /**
   * \brief Provides a transformation function for the nth derivative of the
   * thermal contraction t = t * d^n/d(-beta)^n 1/(1+exp(-/+eps*beta)).
   * For -/+ is used for particles/holes, respectively.
   **/
  template <typename F=real>
  class ThermalContraction {
  public:
    ThermalContraction(
      const real beta_, const bool particle, const int n = 0
    ): beta(beta_), sign(particle ? -1.0 : +1.0), a(n) {
      if (n > 0) {
        std::vector<int64_t> nextA(n);
        a[0] = 1;
        for (int m(3); m <= n+1; ++m) {
          nextA[0] = 1;
          for (int k(1); k < m-1; ++k) {
            // build coefficients for nth derivative
            nextA[k] = (k+1)*a[k] - (m-k-1)*a[k-1];
          }
          a = nextA;
        }
      }
    }
    void operator()(const real eps, F &t) {
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
        t *= std::pow(sign*eps,a.size()) * y;
      }
    }
  protected:
    real beta, sign;
    std::vector<int64_t> a;
  };
}

#endif

