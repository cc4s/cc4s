/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef THERMAL_CLUSTER_DOUBLES_ALGORITHM_DEFINED 
#define THERMAL_CLUSTER_DOUBLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <string>
#include <ctf.hpp>
#include <tcc/DryTensor.hpp>

namespace cc4s {
  /**
   * \brief Provides the functionality for a finite temperature algorithm with
   * only doubles amplitudes. It updates the correlation free energy from
   * the amplitudes.
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

    /**
     * \brief Defines the default number of samples in imaginary time, 32.
     */
    static int constexpr DEFAULT_SAMPLES = 32;

  protected:
    /**
     * \brief The current correlation energy from all diagrams that were
     * closed up to the current moment in imaginart time.
     **/
    CTF::Scalar<> *directEnergy, *exchangeEnergy;

    /**
     * \brief The thermal amplitudes of doubles diagrams of the
     * current and the subsequent moment in imaginary time.
     */
    CTF::Tensor<> *Tabij[2];

    /**
     * \brief The eigenenergy difference of between the state a and i.
     * This is used to propagate all states in imaginary time.
     **/
    CTF::Tensor<> *Dai;

    /**
     * \brief Inverse temperature \f$\beta=1/k_{\rm B}T\f$, where
     * \f$k_{\rm B}T\f$ is given in the same unit as the eigenenergies
     * \f$\varepsilon_p\f$.
     **/
    double beta;

    /**
     * \brief The number of imaginary frequency intervals between \f$0\f$ and
     * \f$\beta\f$ for the approximate solution of the amplitude integral
     * equation.
     **/
    int samples;

    /**
     * \brief Advances the correlation energy and the amplitudes from the
     * the current moment in imaginary time to the next sample according
     * to the concrete implementation.
     * \param[in] i sample number
     */
    virtual void update(int n) = 0;

    /**
     * \brief Performs a dry run of an update according to the concrete
     * algorithm.
     * The base class does not perform accounting and writes a warning about
     * that.
     */
    virtual void dryUpdate();
  };


  class ImaginaryTimePropagation {
  protected:
    ImaginaryTimePropagation(
      double DTau_
    ): DTau(DTau_) {
    }
    double DTau;
  };

  /**
   * \brief Offers a transform method for the propagation of
   * one particle hole pair without interaction
   * within the imgainary time interval of length DTau.
   **/
  class FreePHImaginaryTimePropagation: public ImaginaryTimePropagation {
  public:
    FreePHImaginaryTimePropagation(
      double DTau_
    ): ImaginaryTimePropagation(DTau_) {
    }
    /**
     * \brief Transforms the amplitude P given the energy differences
     * of the particle hole pair.
     * \param[in] Dai The energy difference of the pair propagating.
     * \param[inout] P The amplitude to transform according to the propagation.
     **/
    void operator ()(double Dai, double &P) {
      P *= std::exp(-Dai*DTau);
    }
  };

  /**
   * \brief Offers a transform method for the propagation of
   * two particle hole pairs without interaction
   * within the imgainary time interval of length DTau.
   **/
  class FreePPHHImaginaryTimePropagation: public ImaginaryTimePropagation {
  public:
    FreePPHHImaginaryTimePropagation(
      double DTau_
    ): ImaginaryTimePropagation(DTau_) {
    }
    /**
     * \brief Transforms the amplitude P given the energy differences
     * of the particle hole pairs.
     * \param[in] Dai The energy difference of the first pair propagating.
     * \param[in] Dbj The energy difference of the second pair propagating.
     * \param[inout] P The amplitude to transform according to the propagation.
     **/
    void operator ()(double Dai, double Dbj, double &P) {
      P *= std::exp(-(Dai+Dbj)*DTau);
    }
  };

  /**
   * \brief Offers a transform method for the propagation of
   * two particle hole pairs from or to the Coulomb interaction
   * at tau within the imgainary time interval of length DTau.
   **/
  class PPHHImaginaryTimePropagation: public ImaginaryTimePropagation {
  public:
    PPHHImaginaryTimePropagation(
      double DTau_
    ): ImaginaryTimePropagation(DTau_) {
    }
    /**
     * \brief Transforms the amplitude P given the energy differences
     * of two particle hole pairs, which are either both incoming or both
     * outgoing.
     * \param[in] Dai The energy difference of the first pair propagating.
     * \param[in] Dbj The energy difference of the second propagating.
     * \param[inout] P The amplitude to transform according to the propagation.
     **/
    void operator ()(double Dai, double Dbj, double &P) {
      const double Delta(Dai+Dbj);
      const double DeltaDTau(Delta*DTau);
      // use the first order approximation around Delta=zero if applicable
      // to avoid division by small numbers
      P *= (DeltaDTau*DeltaDTau * DTau > 6e-15) ?
        (1 - std::exp(-DeltaDTau)) / Delta :
        DTau * (1 - DeltaDTau*DTau/2);
    }
  };

  /**
   * \brief Offers a transform method for the propagation of
   * one particle hole pair to and one pair from the Coulomb interaction
   * at tau within the imgainary time interval of length DTau.
   **/
  class HPPHImaginaryTimePropagation: public ImaginaryTimePropagation {
  public:
    HPPHImaginaryTimePropagation(
      double DTau_
    ): ImaginaryTimePropagation(DTau_) {
    }
    /**
     * \brief Transforms the amplitude P given the energy differences
     * of the in- and outgoing particle hole pairs.
     * \param[in] Dai The energy difference of the pair propagating to the
     *                interaction.
     * \param[in] Dbj The energy difference of the pair propagating from the
     *                interaction.
     * \param[inout] P The amplitude to transform according to the propagation.
     **/
    void operator ()(double Dai, double Dbj, double &P) {
      const double meanD((Dbj+Dai)*0.5);
      const double Delta((Dbj-Dai)*0.5);
      const double meanP(std::exp(-meanD*DTau));
      const double DeltaDTau(Delta*DTau);
      const double DeltaDTauSquared(DeltaDTau*DeltaDTau);
      // use the second order approximation around Delta=zero if applicable
      // to avoid division by small numbers
      P *= meanP * (
        (DeltaDTauSquared*DeltaDTauSquared * DTau > 1.20e-13) ?
          std::sinh(DeltaDTau) / Delta :
          DTau * (1 + DeltaDTauSquared*DTau/6)
      );
    }
  };
}

#endif

