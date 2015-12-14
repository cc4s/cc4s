/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COULOMB_RPA_DEFINED 
#define COULOMB_RPA_DEFINED

#include <Algorithm.hpp>
#include <Options.hpp>

namespace cc4s {
  class CoulombRpa: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CoulombRpa);
    CoulombRpa(
      std::vector<Argument> const &argumentList
    );
    virtual ~CoulombRpa();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new CoulombRpa(argumentList);
    }

  protected:
    void iterateCoulombRpa();

   // static properties, accessible from everywhere
    static Options *options;
  };
}

#endif

