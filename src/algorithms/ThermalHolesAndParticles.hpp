/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef THERMAL_HOLES_AND_PARTICLES_DEFINED
#define THERMAL_HOLES_AND_PARTICLES_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  /**
   * \brief Calculates the chemical potential \f$\mu\f$ such that the
   * thermal expectation value of the electron number operator
   * \f$\langle \sum_p\hat p^\dagger \hat p \rangle_0\f$
   * in the reference system equals the given number of electrons.
   * The thermal particle- and hole eigenenergies are then returned
   * relative to \f$\mu\f$.
   * Additionally, the particle- and hole occupancies are calculated
   * given by \f$f^+_p = \langle \hat p \hat p\dagger \rangle_0\f$ and
   * \f$f^-_p = \langle \hat p^\dagger \hat p \rangle_0\f$, respectively.
   * The index set of particles and holes is truncated where the
   * occupancies lead to negligible contributions according to
   * \f$\varepsilon_p-\mu \lessgtr \pm
   *   (32k_{\rm B}T(\varepsilon_{\rm max}-\varepsilon_{\rm min}))^{1/2} ) \f$.
   */
  class ThermalHolesAndParticles: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ThermalHolesAndParticles);
    ThermalHolesAndParticles(
      std::vector<Argument> const &argumentList
    );
    virtual ~ThermalHolesAndParticles();
    virtual void run();
    virtual void dryRun();
  protected:
    void orderStates();
    void determineChemicalPotential();
    void determineNumberOfElectrons();
    void defineThermalHolesAndParticles();
    void determineThermalOccupancies();
    void computeParticleHoleOverlap();

    CTF::Vector<> *epsp;
    std::vector<std::pair<double, int>> eigenStates;
    double kT, mu;
  };
}

#endif

