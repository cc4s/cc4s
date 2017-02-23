/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRCCD_DENSITIES_DEFINED 
#define DRCCD_DENSITIES_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  /**
   * \brief Implements the iteration routine for the Drccd method. Calculates the
   * amplitudes \f$T_{ab}^{ij}\f$ from the Coulomb Integrals \f$V_{ij}^{ab}\f$
   * in a \f$ \mathcal{O}(N^{6}) \f$ implementation.
   */
  class DrccdDensities: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DrccdDensities);
    DrccdDensities(
      std::vector<Argument> const &argumentList
    );
    virtual ~DrccdDensities();

    virtual void run();
    virtual void dryRun();

  protected:
    template<typename T, typename MT>
    void run(T *ctfDabij, const bool dry = false);

    static constexpr int DEFAULT_MAX_ITERATIONS = 16;
  };
}

#endif

