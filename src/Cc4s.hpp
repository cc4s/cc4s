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
    Cc4s(CTF::World *world, Options *options);
    ~Cc4s();
    void run();
    void testSymmetries();

  protected:
    void iterateMp2();
    void iterateCcsd();
    void iterateRpa();
    void iterateRccd();
    void iterateRccsd();

    /**
     * \brief Read Fourier transformed overlap densities and eigenergies from
     * disk and calculate all necessary quantities.
     */
    void readFTOD();
    /**
     * \deprecated
     */
    void add_Vxyef_T21efij(CTF::Tensor<> &Zabij, CTF::Tensor<> &T21);

    CTF::World *world;
    bool profile;

    Options *options;
    Chi *chiReal, *chiImag;
    CoulombIntegrals *V;
    Amplitudes *T;
};

#endif

