/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CHI_DEFINED
#define CHI_DEFINED

#include "PerturbationTensor.hpp"
#include "Options.hpp"
#include <ctf.hpp>

class Chi: public PerturbationTensor {
  public:
    Chi(CTF::World *world, Options const &options, int seed=0);
    virtual ~Chi();

    virtual CTF::Idx_Tensor get(char const *stdIndexMap, char const *indexMap);

    CTF::Tensor<> getSlice(int a);
    void readRandom(CTF::Tensor<> *tensor, int seed);

    CTF::Tensor<> *ab, *ai, *ij;

    // read from disk
    void read();
};

#endif

