/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CC4S_DEFINED
#define CC4S_DEFINED

#include <ctf.hpp>
#include "Options.hpp"
#include "Chi.hpp"
#include "CoulombIntegrals.hpp"
#include "Amplitudes.hpp"

class Cc4s {
  public:
    Cc4s(CTF::World *world, Options const &options);
    ~Cc4s();
    void run();

  protected:
    void iterateAmplitudes();
    void add_Vxyef_T21efij(CTF::Tensor<> &Zabij, CTF::Tensor<> &T21);

    CTF::World *world;
    int no, nv, nG, niter;
    bool profile;

    Chi *chi;
    CoulombIntegrals *V;
    Amplitudes *T;
};

#endif

