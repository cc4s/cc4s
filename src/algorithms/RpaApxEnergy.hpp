/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef RPA_APX_ENERGY_DEFINED 
#define RPA_APX_ENERGY_DEFINED

#include <algorithms/Algorithm.hpp>

#include <util/SharedPointer.hpp>

namespace cc4s {
  class RpaApxEnergy: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(RpaApxEnergy);
    RpaApxEnergy(
      std::vector<Argument> const &argumentList
    );
    virtual ~RpaApxEnergy();

    virtual void run();
  protected:
    void computeFrequencyGrid();
  };
}

#endif

