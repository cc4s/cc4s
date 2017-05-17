/*Copyright (c) 2017, Andreas Grueneis, Felix Hummel and Alejandro Gallo, all rights reserved.*/
#ifndef DRCCD_DENSITIES_DEFINED
#define DRCCD_DENSITIES_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  /**
   * \brief Implements the iteration routine for the Drccd method. Calculates the
   * amplitudes \f$T_{ab}^{ij}\f$ from the Coulomb Integrals \f$V_{ij}^{ab}\f$
   * in a \f$ \mathcal{O}(N^{6}) \f$ implementation.
   */
  class DrccdEquationOfMotion: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DrccdEquationOfMotion);
    DrccdEquationOfMotion(
      std::vector<Argument> const &argumentList
    );
    virtual ~DrccdEquationOfMotion();

    virtual void run();

  protected:
    static constexpr int DEFAULT_MAX_ITERATIONS = 16;
    double energyShift;

    void determineEnergyShift();
  };
}

#endif

