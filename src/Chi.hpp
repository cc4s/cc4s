/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CHI_DEFINED
#define CHI_DEFINED

#include "Options.hpp"
#include <ctf.hpp>

enum ChiPart {
  GAB, GAI, GIJ
};

class Chi {
  public:
    Chi(CTF::World *world, Options const &options);
    ~Chi();
    CTF::Tensor<> &get(ChiPart part);
    CTF::Tensor<> getSlice(ChiPart, int a);

    void readRandom(CTF::Tensor<> *tensor, int seed);
  private:
    CTF::Tensor<> *ab, *ai, *ij;

    // read from disk
    void read();
};

#endif

