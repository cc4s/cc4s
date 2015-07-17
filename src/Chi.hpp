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

  private:
    CTF::Tensor<> *ab, *ai, *ij;
    void readRandom(CTF::Tensor<> *tensor);

    // read from disk
    void read();
};

#endif

