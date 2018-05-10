/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef THERMAL_CLUSTER_DOUBLES_ALGORITHM_DEFINED 
#define THERMAL_CLUSTER_DOUBLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <util/SharedPointer.hpp>
#include <string>
#include <ctf.hpp>
#include <tcc/DryTensor.hpp>

namespace cc4s {
  /**
   * \brief Provides the functionality for a finite temperature algorithm with
   * only doubles amplitudes.
   */
  class ThermalClusterDoublesAlgorithm: public Algorithm {
  public:
    ThermalClusterDoublesAlgorithm(
      std::vector<Argument> const &argumentList
    );
    virtual ~ThermalClusterDoublesAlgorithm();
    /**
     * \brief Calculates the energy of this ThermalClusterDoubles algorithm
     */
    virtual void run();

    /**
     * \brief Performs a Dry Run
     */
    virtual void dryRun();
    /**
     * \brief Returns the abbreviation of the concrete algorithm, without
     * the leading "Thermal", e.g. "Ccd" for "ThermalCcdEnergy"
     */
    virtual std::string getAbbreviation() = 0;

    static constexpr int DEFAULT_MAX_ITERATIONS = 12;

  protected:
    /**
     * \brief doubles amplitudes on the imaginary time grid
     **/
    std::vector<PTR(CTF::Tensor<real>)> Tabijn;

    /**
     * \brief The eigenenergy difference of between the state a,b and i,j.
     * This is used to propagate all states in imaginary time.
     **/
    PTR(CTF::Tensor<>) Dabij;

    /**
     * \brief Inverse temperature \f$\beta=1/k_{\rm B}T\f$, where
     * \f$k_{\rm B}T\f$ is given in the same unit as the eigenenergies
     * \f$\varepsilon_p\f$.
     **/
    real beta;


    /**
     */
    virtual void applyHamiltonian(
      const CTF::Tensor<real> &T0abij,
      CTF::Tensor<real> &T1abij,
      real DTau
    ) = 0;

    std::string getCapitalizedAbbreviation();
    std::string getAmplitudeIndices(CTF::Tensor<> &T);
    void fetchDelta(CTF::Tensor<> &Delta);
    void thermalContraction(CTF::Tensor<> &T);

    class ImaginaryTimeTransform {
    protected:
      ImaginaryTimeTransform(
        real DTau_
      ): DTau(DTau_) {
      }
      real DTau;
    };

    class FreeImaginaryTimePropagation: public ImaginaryTimeTransform {
    public:
      FreeImaginaryTimePropagation(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      void operator ()(const real Delta, real &T) const {
        T *= std::exp(-Delta*DTau);
      }
    };

    // convolves the T^I(tau_m-1) contribution to T'^I(tau_m)
    // = f^(J\I) * (1-(1+DeltaJ*Tau)*exp(-DeltaJ*Tau))/(DeltaJ^2*Tau)
    class BeginConvolution: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-6;
      BeginConvolution(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      // \brief Requires P = log(f^(J\I))
      void operator ()(const real DeltaJ, real &P) {
        const real x(DTau*DeltaJ);
        if (x > SMALL) {
          P += std::log1p(-(1+x)*std::exp(-x)) - std::log(+x);
          P = std::exp(P) / (+DeltaJ);
        } else if (x < -SMALL) {
          P += std::log1p(-(1+x)*std::exp(-x)) - std::log(-x);
          P = std::exp(P) / (-DeltaJ);
        } else {
          P = std::exp(P) * DTau * ( 0.5 - x*(1./3 + x*(1./8 - x/30)) );
        }
      }
    };
  };

    // convolves the T^I(tau_m) contribution to T'^I(tau_m)
    // = f^(J/I) * (DeltaJ*Tau-(1-exp(-DeltaJ*Tau)))/(DeltaJ^2*Tau)
    class EndConvolution: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-6;
      static constexpr real LARGE = 42.0;
      EndConvolution(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      // \brief Requires P = log(f^(J\I))
      void operator ()(const real DeltaJ, real &P) {
        const real x(DTau*DeltaJ);
        if (std::abs(x) < SMALL) {
          P = std::exp(P) * DTau * ( 0.5 - x*(1./6 + x*(1./24 - x/120)) );
        } else if (x > -LARGE) {
          P = std::exp(P) * (1-(1-std::exp(-x))/x) / DeltaJ;
        } else
          P += -x - std::log(-x);
          P = std::exp(P) / (-DeltaJ);
        }
      }
    };
  };
}

#endif

