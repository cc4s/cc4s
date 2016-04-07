/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED 
#define CCD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/ClusterDoublesAlgorithm.hpp>

namespace cc4s {
  // this algorithm is now based on the ClusterDoublesAlgorithm
  // inheriting its iteration and slicing functionality.
  // Only the abstract (left out) methods getAbbreviation and iterate have
  // to be implemented.
  class CcdEnergyFromCoulombIntegrals: public ClusterDoublesAlgorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcdEnergyFromCoulombIntegrals);
    CcdEnergyFromCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~CcdEnergyFromCoulombIntegrals();
    virtual std::string getAbbreviation() { return "Ccd"; }

  protected:
    virtual void iterate(int i);
    virtual void iterateBartlett(int i);
  };
}

#endif

