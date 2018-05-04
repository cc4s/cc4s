/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef THERMAL_DIRECTOR_RING_COUPLED_CLUSTER_DOUBLES_DEFINED 
#define THERMAL_DIRECTOR_RING_COUPLED_CLUSTER_DOUBLES_DEFINED

#include <algorithms/ThermalClusterDoublesAlgorithm.hpp>

namespace cc4s {
  /**
   * \brief Implements an iterative solution of the Direct Rin
   * Coupled Cluster amplitude integral equation at finite temperatures.
   */
  class ThermalDirectRingCoupledClusterDoubles:
    public ThermalClusterDoublesAlgorithm
  {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ThermalDirectRingCoupledClusterDoubles);
    ThermalDirectRingCoupledClusterDoubles(
      std::vector<Argument> const &argumentList
    );
    virtual ~ThermalDirectRingCoupledClusterDoubles();

    virtual std::string getAbbreviation() {
      return "Drccd";
    }

  protected:
    virtual void applyHamiltonian(
      const CTF::Tensor<real> &T0abij,
      CTF::Tensor<real> &T1abij,
      real DTau
    );
  };
}

#endif

