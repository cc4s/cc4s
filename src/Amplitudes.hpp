/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef AMPLITUDES_DEFINED
#define AMPLITUDES_DEFINED

#include <ctf.hpp>
#include "Exception.hpp"
#include "CoulombIntegrals.hpp"

class Amplitudes {
  public:
    Amplitudes(CoulombIntegrals *V);
    CTF::Tensor<> &get(Part part);

  protected:
    CTF::Tensor<> *ai, *abij;
};

#endif

