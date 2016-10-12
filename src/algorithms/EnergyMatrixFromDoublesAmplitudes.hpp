/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef ENERGY_MATRIX_FROM_DOUBLES_AMPLITUDES_DEFINED
#define ENERGY_MATRIX_FROM_DOUBLES_AMPLITUDES_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class EnergyMatrixFromDoublesAmplitudes: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(EnergyMatrixFromDoublesAmplitudes);
    EnergyMatrixFromDoublesAmplitudes(std::vector<Argument> const &argumentList);
    virtual ~EnergyMatrixFromDoublesAmplitudes();
    /**
     * \brief calculates the energy matrix \f$E_G^{G'}\f$ from the given
     * doubles amplitudes \f$T_{ij}^{ab}\f$ such that the direct term
     * of the energy \f$T_{ij}^{ab}V_{ab}^{ij}\f$ is given by the trace
     * of the energy matrix \f$E_G^G\f$.
     * This algorithm is intended for testing purposes only, since this
     * calculation is usually done in the DFT/HF electronic structure code
     * in order to reduces the size of the Coulomb vertex.
     */
    virtual void run();
    /**
     * \brief Dry run for calculating the energy matrix from the given
     * doubles amplitdues.
     */
    virtual void dryRun();
  };
}

#endif

