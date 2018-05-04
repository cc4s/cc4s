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
  };

/*
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

  class ImaginaryTimeConvolution: public ImaginaryTimeTransform {
  public:
    static constexpr real EPSILON = 1e-6;
    ImaginaryTimeConvolution(
      real DTau_
    ): ImaginaryTimeTransform(DTau_) {
    }
     // \brief Requires P = -DeltaJ*DTau + log(f^J) + log(f^(K\J)) - log(f^I)
    void operator ()(const real DeltaK, real &P) {
      if (DeltaK*DTau > EPSILON) {
        P -= std::log1p(-std::exp(-DeltaK*DTau));
        P = std::exp(P) / (+DeltaK);
      } else if ((DeltaK*DTau < EPSILON) {
        P -= std::log1p(-std::exp(+DeltaK*DTau)) + DeltaK*DTau;
        P = std::exp(P) / (-DeltaK);
      } else {
        P -= DeltaK*DTau/2 + (DeltaK*DTau)*(DeltaK*DTau)/24;
        P = std::exp(P) * DTau;
      }
    }
  };
*/
}

#endif

