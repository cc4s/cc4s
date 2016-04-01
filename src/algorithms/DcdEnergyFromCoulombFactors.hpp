/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DCD_ENERGY_FROM_COULOMB_FACTORS_DEFINED 
#define DCD_ENERGY_FROM_COULOMB_FACTORS_DEFINED

#include <algorithms/Algorithm.hpp>
#include <mixers/Mixer.hpp>

namespace cc4s {
  class DcdEnergyFromCoulombFactors: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DcdEnergyFromCoulombFactors);
    DcdEnergyFromCoulombFactors(
      std::vector<Argument> const &argumentList
    );
    virtual ~DcdEnergyFromCoulombFactors();
    virtual void run();

    static Algorithm *create(std::vector<Argument> const &argumentList) {
      return new DcdEnergyFromCoulombFactors(argumentList);
    }

    static int64_t constexpr DEFAULT_MAX_ITERATIONS = 16;

  protected:
    Mixer<double> *TabijMixer;
    void iterate(int i);
  };
}

#endif

