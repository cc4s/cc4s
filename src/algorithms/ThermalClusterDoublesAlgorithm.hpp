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
    std::vector<PTR(CTF::Tensor<complex>)> Tabijn;
    /**
     * \brief doubles amplitudes on the imaginary frequency grid
     **/
    std::vector<PTR(CTF::Tensor<complex>)> Tabijv;

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
    double beta;


    /**
     * \brief Calculates the residuum for the given imaginary time amplitudes.
     */
    virtual void getResiduum(CTF::Tensor<complex> &Tabij) = 0;

    /**
     * \brief Performs a dry run of an iterate according to the concrete
     * algorithm.
     * The base class does not perform accounting and writes a warning about
     * that.
     */
    virtual void dryIterate();

    std::string getCapitalizedAbbreviation();
    std::string getAmplitudeIndices(CTF::Tensor<> &T);
    void fetchDelta(CTF::Tensor<> &Delta);
    void thermalContraction(CTF::Tensor<> &T);
  };
}

#endif

