/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef RANDOM_ORBITALS_DEFINED
#define RANDOM_ORBITALS_DEFINED

#include <algorithms/Algorithm.hpp>
#include <ctf.hpp>

namespace cc4s {
  class RandomOrbitals: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(RandomOrbitals);
    RandomOrbitals(std::vector<Argument> const &argumentList);
    virtual ~RandomOrbitals();
    virtual void run();
    virtual void dryRun();
  };
}

#endif

