/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef THERMAL_CLUSTER_DOUBLES_ALGORITHM_DEFINED 
#define THERMAL_CLUSTER_DOUBLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <string>
#include <ctf.hpp>
#include <util/DryTensor.hpp>

namespace cc4s {
  /**
   * \brief Provides the functionality for a finite temperature algorithm with
   * only doubles amplitudes. It updates the correlation free energy from
   * the amplitudes.
   */
  class ThermalClusterDoublesAlgorithm: public Algorithm {
  public:
    ClusterDoublesAlgorithm(
      std::vector<Argument> const &argumentList
    );
    virtual ~ClusterDoublesAlgorithm();
    /**
     * \brief Calculates the energy of this ThermalClusterDoubles algorithm
     */
    virtual void run();

    /**
     * \brief Performs a Dry Run
     */
    virtual void dryRun();
    /**
     * \brief Returns the abbreviation of the concrete algorithm, e.g.
     * "ThermalCcd".
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
    Scalar<> *energy;

    /**
     * \brief The current thermal amplitudes of doubles diagrams up to the
     * current moment in imaginary time.
     */
    Tensor<> *Tabij;

    /**
     * \brief Advances the correlation energy and the amplitudes from the
     * the current moment in imaginary time to the next sample according
     * to the concrete implementation.
     * \param[in] i sample number
     */
    virtual void update(int i) = 0;

    /**
     * \brief Performs a dry run of an update according to the concrete
     * algorithm.
     * The base class does not perform accounting and writes a warning about
     * that.
     */
    virtual void dryUpdate();
  };
}

#endif

