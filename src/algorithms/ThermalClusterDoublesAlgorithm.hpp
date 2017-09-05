/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef THERMAL_CLUSTER_DOUBLES_ALGORITHM_DEFINED 
#define THERMAL_CLUSTER_DOUBLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
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

    /**
     * \brief Defines the default recursion length.
     */
    static int constexpr DEFAULT_RECURSION_LENGTH = 2;
    /**
     * \brief Defines the default minimum number of iterations
     */
    static int constexpr DEFAULT_MIN_ITERATIONS = 8;
    /**
     * \brief Defines the default maximum number of iterations
     */
    static int constexpr DEFAULT_MAX_ITERATIONS = 1024;

  protected:
    /**
     * \brief The free energies of all preceeding time scales.
     **/
    std::vector<CTF::Scalar<> *>energies;

    /**
     * \brief The thermal doubles amplitudes of all preceeding time scales.
     */
    std::vector<CTF::Tensor<> *>amplitudes;

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
     * \brief The number of imaginary frequency intervals retained.
     **/
    int recursionLength;

    /**
     * \brief The scaling factor of two successive imaginary time scales.
     **/
    double recursionScaling;

    /**
     * \brief Calculates properties at the next larger imaginary time scale
     * from the properties of the preceding smaller imaginary time scales
     * \param[in] n The decreasing level of energy. At 0 the calcuation is done.
     */
    virtual void iterate(int n) = 0;

    /**
     * \brief Performs a dry run of an iterate according to the concrete
     * algorithm.
     * The base class does not perform accounting and writes a warning about
     * that.
     */
    virtual void dryIterate();

    void initializeRecursion(const int N);
    double recurse(const int n);

    std::string getCapitalizedAbbreviation();
    double getRecursionScaling(const int M);
    double getEnergyScale();
    std::string getAmplitudeIndices(CTF::Tensor<> &T);
    void fetchDelta(CTF::Tensor<> &Delta);
    void thermalContraction(CTF::Tensor<> &T);
  };


  class ImaginaryTimeTransform {
  protected:
    ImaginaryTimeTransform(
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
  class FreePHImaginaryTimeTransform: public ImaginaryTimeTransform {
  public:
    FreePHImaginaryTimeTransform(
      double DTau_
    ): ImaginaryTimeTransform(DTau_) {
    }
    /**
     * \brief Transforms the amplitude P given the energy differences
     * of the particle hole pair.
     * \param[in] Dai The energy difference of the pair propagating.
     * \param[inout] P The amplitude to transform according to the propagation.
     **/
    void operator ()(double Dai, double &P) const {
      P *= std::exp(-Dai*DTau);
    }
  };

  class FreeImaginaryTimePropagation: public ImaginaryTimeTransform {
  public:
    FreeImaginaryTimePropagation(
      double DTau_
    ): ImaginaryTimeTransform(DTau_) {
    }
    /**
     * \brief Returns propagator given the energy sum of the states
     * propagating.
     * \param[in] Delta The propagating energies. Particles count positive.
     **/
    double operator ()(const double Delta) const {
      return std::exp(-Delta*DTau);
    }
  };

  /**
   * \brief Offers a transform method for the propagation of
   * two particle hole pairs without interaction
   * within the imgainary time interval of length DTau.
   **/
  class FreePPHHImaginaryTimeTransform: public ImaginaryTimeTransform {
  public:
    FreePPHHImaginaryTimeTransform(
      double DTau_
    ): ImaginaryTimeTransform(DTau_) {
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
  class PPHHImaginaryTimeTransform: public ImaginaryTimeTransform {
  public:
    PPHHImaginaryTimeTransform(
      double DTau_
    ): ImaginaryTimeTransform(DTau_) {
    }
    /**
     * \brief Transforms the amplitude P given the energy differences
     * of two particle hole pairs, which are either both incoming or both
     * outgoing.
     * \param[in] Dai The energy difference of the first pair propagating.
     * \param[in] Dbj The energy difference of the second propagating.
     * \param[inout] P The amplitude to transform according to the propagation.
     **/
    void operator ()(double Dai, double Dbj, double &P) const {
      const double Delta(Dai+Dbj);
      const double DeltaDTau(Delta*DTau);
      // use the first order approximation around Delta=zero if applicable
      // to avoid division by small numbers
      P *= (DeltaDTau*DeltaDTau * DTau > 6e-15) ?
        (1 - std::exp(-DeltaDTau)) / Delta :
        DTau * (1 - 0.5*DeltaDTau*DTau);
    }
  };

  class SameSideConnectedImaginaryTimePropagation:
    public ImaginaryTimeTransform
  {
  public:
    SameSideConnectedImaginaryTimePropagation(
      double DTau_
    ): ImaginaryTimeTransform(DTau_) {
    }
    /**
     * \brief Returns the propgator for states connected to the interaction
     * given the energy sum of the states propagating. The states must be
     * connected either all from below or all to above.
     * Particles count positive, holes count negative.
     * \param[in] Delta The energy sum of the states propagating.
     **/
    double operator ()(const double Delta) const {
      const double DeltaDTau(Delta*DTau);
      // use the first order approximation around Delta=zero if applicable
      // to avoid division by small numbers
      return (DeltaDTau*DeltaDTau * DTau > 6e-15) ?
        (1 - std::exp(-DeltaDTau)) / Delta :
        DTau * (1 - 0.5*DeltaDTau);
    }
  };

  /**
   * \brief Offers a transform method for the propagation of
   * one particle hole pair to and one pair from the Coulomb interaction
   * at tau within the imgainary time interval of length DTau.
   **/
  class HPPHImaginaryTimeTransform: public ImaginaryTimeTransform {
  public:
    HPPHImaginaryTimeTransform(
      double DTau_
    ): ImaginaryTimeTransform(DTau_) {
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
    void operator ()(double Dai, double Dbj, double &P) const {
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

  class UpDownConnectedImaginaryTimePropagation:
    public ImaginaryTimeTransform
  {
  public:
    UpDownConnectedImaginaryTimePropagation(
      double DTau_
    ): ImaginaryTimeTransform(DTau_) {
    }
    /**
     * \brief Returns the propagator given the energy sums of the
     * states connected from below as well as the energy sums of the states
     * connected to above.
     * Particles count positive, holes count negative.
     * \param[in] D0 The energy sums of the states connected from below.
     * \param[in] D1 The energy sums of the states connected to above.
     * Note that the propagation is symmetric under exchange of D0 and D1.
     **/
    double operator ()(const double D0, const double D1) const {
      const double meanD((D1+D0)*0.5);
      const double Delta((D1-D0)*0.5);
      const double meanP(std::exp(-meanD*DTau));
      const double DeltaDTau(Delta*DTau);
      const double DeltaDTauSquared(DeltaDTau*DeltaDTau);
      // use the second order approximation around Delta=zero if applicable
      // to avoid division by small numbers
      return meanP * (
        (DeltaDTauSquared*DeltaDTauSquared * DTau > 1.20e-13) ?
          std::sinh(DeltaDTau) / Delta :
          DTau * (1 + DeltaDTauSquared*DTau/6)
      );
    }
  };

  class Mp2ImaginaryTimePropagation: public ImaginaryTimeTransform {
  public:
    Mp2ImaginaryTimePropagation(
      const double DTau_
    ): ImaginaryTimeTransform(DTau_) {
    }
    double operator ()(const double Delta) const {
      // use the first order approximation around Delta=zero if applicable
      // to avoid division by small numbers
      const double DeltaDTau( Delta*DTau );
      const double DTauDTau( DTau*DTau );
      const double DeltaDelta( Delta*Delta );
      return DeltaDTau*DeltaDTau*DTauDTau > 2.4e-14 ?
        (std::exp(-DeltaDTau) - 1.0 + DeltaDTau) / DeltaDelta :
        DTau*(0.5*DTau + DeltaDTau*DTau/6);
    }
  };
}

#endif

