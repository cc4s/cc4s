/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DRCCD_ENERGY_FROM_COULOMB_VERTEX_DEFINED 
#define DRCCD_ENERGY_FROM_COULOMB_VERTEX_DEFINED

#include <Algorithm.hpp>
#include <Options.hpp>
#include <util/Complex.hpp>

namespace cc4s {
  class DrccdEnergyFromCoulombVertex: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DrccdEnergyFromCoulombVertex);
    DrccdEnergyFromCoulombVertex(
      std::vector<Argument> const &argumentList
    );
    virtual ~DrccdEnergyFromCoulombVertex();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new DrccdEnergyFromCoulombVertex(argumentList);
    }

  protected:
    void iterate();

    CTF::Tensor<complex> *vabij;
    CTF::Tensor<> *realGammaGai, *imagGammaGai;
  };
}

#endif

