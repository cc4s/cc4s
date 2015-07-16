/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CHI_DEFINED
#define CHI_DEFINED

#include <ctf.hpp>

class Chi {
  public:
    CTF::Tensor<> *ab, *ai, *ij;

    Chi();
    ~Chi();
    // read from disk
    void read();
};

#endif

