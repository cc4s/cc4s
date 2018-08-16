/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
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
    static constexpr int DEFAULT_SINGLES_ENERGY = 1;

  protected:
    /**
     * \brief singles and doubles amplitudes on the imaginary time grid
     **/
    std::vector<PTR(CTF::Tensor<real>)> TFn, TFGn;

    CTF::Tensor<real> *lambdaF;

    /**
     * \brief lambdaF + lambdaG for doubles propagation
     **/
    PTR(CTF::Tensor<real>) lambdaFG;

    /**
     * \brief direct & exchange Vabij in left/right singles eigensystem F/G
     **/
    PTR(CTF::Tensor<real>) VdFG, VxFG;

    /**
     * \brief H_0^a_i in singles eigensystem F
     **/
    PTR(CTF::Tensor<real>) H0F;

    /**
     * \brief sqrt(occupancies)
     **/
    PTR(CTF::Tensor<real>) gi, ga;

    /**
     * \brief Inverse temperature \f$\beta=1/k_{\rm B}T\f$, where
     * \f$k_{\rm B}T\f$ is given in the same unit as the eigenenergies
     * \f$\varepsilon_p\f$.
     **/
    real beta;

    /**
     * \brief eigenvectors of singles part of the Hamiltonian
     **/
    PTR(CTF::Tensor<real>) UaiF;

    /**
     */
    virtual void applyHamiltonian(
      CTF::Tensor<real> &T0F,
      CTF::Tensor<real> &T0FG,
      CTF::Tensor<real> &T1F,
      CTF::Tensor<real> &T1FG,
      const real DTau,
      CTF::Tensor<real> &S1F,
      CTF::Tensor<real> &S1FG
    ) = 0;

    std::string getCapitalizedAbbreviation();
    std::string getAmplitudeIndices(CTF::Tensor<real> &T);
    real getTammDancoffEnergy();
    void computeEnergyContribution(
      CTF::Tensor<real> &T1F, CTF::Tensor<real> &T2FG, const real DTau,
      real &direct, real &exchange, real &singles
    );
    void computeSqrtOccupancies();
    void diagonalizeSinglesHamiltonian();
    void diagonalizeDoublesAmplitudes();

    class ImaginaryTimeTransform {
    protected:
      ImaginaryTimeTransform(
        real DTau_
      ): DTau(DTau_) {
      }
      real DTau;
    };

    class SecondOrderIntegral: public ImaginaryTimeTransform {
    public:
      SecondOrderIntegral(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      void operator ()(const real lambda, real &hh) const {
        const real x(lambda * DTau);
        if (std::abs(x) > 0.25) {
          hh *= DTau * (std::exp(-x) - 1.0 + x) / (x*x);
        } else {
          hh *= DTau/2*(
            1 - x/3*(
              1 - x/4*(
                1 - x/5*(
                  1 - x/6*(
                    1 - x/7*(
                      1 - x/8*(
                        1 - x/9*(
                          1 - x/10*(
                            1 - x/11*(
                              1 - x/12
                            )
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          );
        }
      }
    };

    class ImaginaryTimePropagation: public ImaginaryTimeTransform {
    public:
      ImaginaryTimePropagation(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      void operator ()(const real lambdaFG, real &T) const {
        T *= std::exp(-lambdaFG*DTau);
      }
    };

    // constant V^J=C contribution to T'^J(tau_m),
    // independent of the amplitudes T^I(tau)
    // = int_0^Tau dtau exp(-lambdaFG*(Tau-tau))
    class ConvolutionC: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-2;
      ConvolutionC(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      // \brief Requires T = f^J
      void operator ()(const real lambdaFG, real &T) const {
        const real x(DTau*lambdaFG);
        if (std::abs(x) < SMALL) {
          T *= DTau * (
            1. - x*(
              1./2 - x*(
                1./6 - x*(1./24 - x*(1./120 - x*(1./720 - x/5040)))
              )
            )
          );
        } else {
          T *= (1-std::exp(-x)) / lambdaFG;
        }
      }
    };

    // linear T^I(tau_m-1)=T0 contribution to T'^J(tau_m)
    // = int_0^Tau dtau (Tau-tau)/Tau * exp(-lambdaFG*(Tau-tau))
    class Convolution0: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-1;
      Convolution0(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      // \brief Requires T = f^(J\I)
      void operator ()(const real lambdaFG, real &T) {
        const real x(DTau*lambdaFG);
        if (std::abs(x) < SMALL) {
          T *= DTau * (
            1./2 - x*(
              1./3 - x*(
                1./8 - x*(
                  1./30 - x*(1./144 - x*(1./840 - x*(1./5760 - x/45360)))
                )
              )
            )
          );
        } else {
          T *= (1-(x+1)*std::exp(-x)) / (x*lambdaFG);
        }
      }
    };

    // linear T^I(tau_m)=T1 contribution to T'^J(tau_m)
    // = int_0^Tau dtau tau/Tau * exp(-lambdaFG*(Tau-tau))
    class Convolution1: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-1;
      Convolution1(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      // \brief Requires T = f^(J\I)
      void operator ()(const real lambdaFG, real &T) {
        const real x(DTau*lambdaFG);
        if (std::abs(x) < SMALL) {
          T *= DTau * (
            1./2 - x*(
              1./6 - x*(
                1./24 - x*(
                  1./120 - x*(1./720 - x*(1./5040 - x*(1./40320 - x/362880)))
                )
              )
            )
          );
        } else {
          T *= (std::exp(-x)+x-1) / (x*lambdaFG);
        }
      }
    };

    // TODO: improve low x behaviour
    // quadratic T^I1(tau_m-1)*T^I2(tau_m-1)=T0*T0 contribution to T'^J(tau_m)
    // = int_0^Tau dtau (Tau-tau)/Tau * (Tau-tau)/Tau * exp(-lambdaFG*(Tau-tau))
    class Convolution00: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-6;
      Convolution00(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      // \brief Requires T = f^(J\I)
      void operator ()(const real lambdaFG, real &T) {
        const real x(DTau*lambdaFG);
        if (std::abs(x) < SMALL) {
          T *= DTau * ( 1./3 - x*(1./4 - x*(1./10 - x/36)) );
        } else {
          T *= (2-((2+2*x+x*x)*std::exp(-x))) / (x*x*lambdaFG);
        }
      }
    };

    // quadratic T^I1(tau_m-1)*T^I2(tau_m)=T0*T1 contribution to T'^J(tau_m)
    // = int_0^Tau dtau (Tau-tau)/Tau * tau/Tau * exp(-lambdaFG*(Tau-tau))
    class Convolution01: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-6;
      Convolution01(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      void operator ()(const real lambdaFG, real &T) {
        const real x(DTau*lambdaFG);
        if (std::abs(x) < SMALL) {
          T *= DTau * ( 1./6 - x*(1./12 - x*(1./40 - x/180)) );
        } else {
          T *= ((x-2)+((x+2)*std::exp(-x))) / (x*x*lambdaFG);
        }
      }
    };

    // quadratic T^I1(tau_m)*T^I2(tau_m)=T1*T1 contribution to T'^J(tau_m)
    // = int_0^Tau dtau tau/Tau * tau/Tau * exp(-lambdaFG*(Tau-tau))
    class Convolution11: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-6;
      Convolution11(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      void operator ()(const real lambdaFG, real &T) {
        const real x(DTau*lambdaFG);
        if (std::abs(x) < SMALL) {
          T *= DTau * ( 1./3 - x*(1./12 - x*(1./60 - x/360)) );
        } else {
          T *= ((2-2*x+x*x)-(2*std::exp(-x)))/ (x*x*lambdaFG);
        }
      }
    };
  };
}

#endif

