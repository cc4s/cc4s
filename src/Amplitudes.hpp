/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef AMPLITUDES_DEFINED
#define AMPLITUDES_DEFINED

#include <PerturbationTensor.hpp>
#include <CoulombIntegrals.hpp>
#include <ctf.hpp>

namespace cc4s {
  class Amplitudes: public PerturbationTensor {
    public:
      Amplitudes(CoulombIntegrals *V);
      virtual ~Amplitudes();
      virtual CTF::Idx_Tensor get(char const *stdIndexMap, char const *indexMap);

      CTF::Tensor<> *ai, *abij;
  };
}

#endif

