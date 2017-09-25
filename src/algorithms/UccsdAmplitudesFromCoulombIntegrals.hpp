/*Copyright (c) 2017, Andreas Grueneis, Felix Hummel and Alejandro Gallo, all rights reserved.*/
#ifndef MP2_EOM_DEFINED
#define MP2_EOM_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  /**
   * \brief Implements the iteration routine for the Mp2 method. Calculates the
   * amplitudes \f$T_{ab}^{ij}\f$ from the Coulomb Integrals \f$V_{ij}^{ab}\f$
   * in a \f$ \mathcal{O}(N^{6}) \f$ implementation.
   */
  class UccsdAmplitudesFromCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(UccsdAmplitudesFromCoulombIntegrals);
    UccsdAmplitudesFromCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~UccsdAmplitudesFromCoulombIntegrals();

    virtual void run();

  protected:
    static constexpr int DEFAULT_MAX_ITERATIONS = 16;

  };
}

#endif

