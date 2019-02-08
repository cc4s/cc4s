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
    static constexpr int DEFAULT_SINGLES_ENERGY = 0;
    static constexpr int DEFAULT_SINGLES = 0;

  protected:
    std::vector<real> taus;

    CTF::Tensor<real> *lambdaF;

    /**
     * \brief lambdaF + lambdaG for doubles propagation
     **/
    PTR(CTF::Tensor<real>) lambdaFG;

    /**
     * \brief direct & exchange Vabij in left/right singles-mode basis F/G
     **/
    PTR(CTF::Tensor<real>) VdFG, VxFG;

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
      CTF::Tensor<real> &T0FG,
      CTF::Tensor<real> &T1FG,
      const real DTau,
      CTF::Tensor<real> &S1FG
    ) = 0;

    std::string getCapitalizedAbbreviation();
    std::string getAmplitudeIndices(CTF::Tensor<real> &T);
    cc4s::real getTammDancoffEnergy();
    cc4s::real getZeroTDrccd();
    void setupImaginaryTimeGrid();
    void iterateAmplitudeSamples();
    void computeEnergyContribution(
      CTF::Tensor<real> &SFG, const real DTau,
      real &direct, real &exchange
    );
    void computeSqrtOccupancies();
    void diagonalizeSinglesHamiltonian();
    void propagateAmplitudes(
      CTF::Tensor<real> &SFG,
      const std::function<void(real, real &)> &propagator
    );
    void diagonalizeDoublesAmplitudes(
      CTF::Tensor<real> &TFG
    );

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
      static constexpr real SMALL = 0.25;
      ConvolutionC(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      // \brief Requires T = f^J
      void operator ()(const real lambdaFG, real &T) const {
        const real x(DTau*lambdaFG);
        if (std::abs(x) < SMALL) {
          T *= DTau * (
            1 - x/2*(
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
      static constexpr real SMALL = 1e-2;
      Convolution00(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      // \brief Requires T = f^(J\I)
      void operator ()(const real lambdaFG, real &T) {
        const real x(DTau*lambdaFG);
        if (std::abs(x) < SMALL) {
          T *= DTau * (
            1./3 - x*(
              1./4 - x*(
                1./10 - x*(
                  1./36 - x*(
                    1./168 - x*(
                      1./960 - x*(
                        1./6480 - x*(
                          1./50400 - x*(
                            1./443520
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          );
        } else {
          T *= (2-((2+2*x+x*x)*std::exp(-x))) / (x*x*lambdaFG);
        }
      }
    };

    // quadratic T^I1(tau_m-1)*T^I2(tau_m)=T0*T1 contribution to T'^J(tau_m)
    // = int_0^Tau dtau (Tau-tau)/Tau * tau/Tau * exp(-lambdaFG*(Tau-tau))
    class Convolution01: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-2;
      Convolution01(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      void operator ()(const real lambdaFG, real &T) {
        const real x(DTau*lambdaFG);
        if (std::abs(x) < SMALL) {
          T *= DTau * (
            1./6 - x*(
              1./12 - x*(
                1./40 - x*(
                  1./180 - x*(
                    1./1008 - x*(
                      1./6720 - x*(
                        1./51840 - x*(
                          1./453600 - x*(
                            1./4435200
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          );
        } else {
          T *= ((x-2)+((x+2)*std::exp(-x))) / (x*x*lambdaFG);
        }
      }
    };

    // quadratic T^I1(tau_m)*T^I2(tau_m)=T1*T1 contribution to T'^J(tau_m)
    // = int_0^Tau dtau tau/Tau * tau/Tau * exp(-lambdaFG*(Tau-tau))
    class Convolution11: public ImaginaryTimeTransform {
    public:
      static constexpr real SMALL = 1e-2;
      Convolution11(
        real DTau_
      ): ImaginaryTimeTransform(DTau_) {
      }
      void operator ()(const real lambdaFG, real &T) {
        const real x(DTau*lambdaFG);
        if (std::abs(x) < SMALL) {
          T *= DTau * (
            1./3 - x*(
              1./12 - x*(
                1./60 - x*(
                  1./360 - x*(
                    1./2520 - x*(
                      1./20160 - x*(
                        1./181440 - x*(
                          1./1814400 - x*(
                            1./19958400
                          )
                        )
                      )
                    )
                  )
                )
              )
            )
          );
        } else {
          T *= ((2-2*x+x*x)-(2*std::exp(-x)))/ (x*x*lambdaFG);
        }
      }
    };
  };
}

#endif

