/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DCSD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED 
#define DCSD_ENERGY_FROM_COULOMB_INTEGRALS_DEFINED

#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>

namespace cc4s {
  class DcsdEnergyFromCoulombIntegrals: public ClusterSinglesDoublesAlgorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DcsdEnergyFromCoulombIntegrals);
    DcsdEnergyFromCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~DcsdEnergyFromCoulombIntegrals();

    virtual std::string getAbbreviation() { return "Dcsd"; }

  protected:
    virtual void iterate(int i);
  };
}

#endif

