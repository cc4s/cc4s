#ifndef MP2_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED 
#define MP2_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  /**
   * \brief Caclulates MP2 energy from the Coulomb Integrals \f$V_{ij}^{ab}.
   */
  class Mp2EnergyFromCoulombIntegrals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(Mp2EnergyFromCoulombIntegrals)
    /**
     * \brief Calculates MP2 energy from Coulomb integrals Vabij
     */
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
  protected:
    template <typename F, typename TE>
    Ptr<MapNode> calculateMp2Energy(const Ptr<MapNode> &arguments);
  };
}

#endif

