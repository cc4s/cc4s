/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CHI_DEFINED
#define CHI_DEFINED

#include "PerturbationTensor.hpp"
#include <ctf.hpp>

class Chi: public PerturbationTensor {
  public:
    Chi(int nG, int no, int nv);
    virtual ~Chi();
    virtual CTF::Idx_Tensor get(char const *stdIndexMap, char const *indexMap);
    CTF::Tensor<> getSlice(int pStart, int pEnd, int qStart, int qEnd);

    int nG, no, nv;
    CTF::Tensor<> *gpq;
};

#endif

