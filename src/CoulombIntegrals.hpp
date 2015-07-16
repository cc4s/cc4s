/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COULOMB_INTEGRALS_DEFINED
#define COULOMB_INTEGRALS_DEFINED

#include <ctf.hpp>

class CoulombIntegrals {
  public:
    CTF::Tensor<> *a, *i, *ai, *abij;
// NOTE: only for testing
    CTF::Tensor<> *abcd;

    CoulombIntegrals();
    void calculate();
    void calculate_xycd(CTF::Tensor<> &xycd, int a, int b);
};

#endif

