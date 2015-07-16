/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "Chi.hpp"
using namespace CTF;

#include "cc4s.hpp"

Chi::Chi() {
  // keep the chi tensors in the memory for now
  {
    int lens[] = {cc4s::nG, cc4s::nv, cc4s::nv};
    int syms[] = {NS, SY, NS};
    ab = new Tensor<>(3, lens, syms, *cc4s::world, "Xab", cc4s::profile);
  }
  {
    int lens[] = {cc4s::nG, cc4s::nv, cc4s::no};
    int syms[] = {NS, NS, NS};
    ai = new Tensor<>(3, lens, syms, *cc4s::world, "Xai", cc4s::profile);
  }
  {
    int lens[] = {cc4s::nG, cc4s::no, cc4s::no};
    int smys[] = {NS, NS, NS};
    ij = new Tensor<>(3, lens, smys, *cc4s::world, "Xij", cc4s::profile);
  }
}

Chi::~Chi() {
  delete ab; delete ai; delete ij;
}

/**
 * \brief Reads the chi amplitudes from disk
 */
void Chi::read() {
  // TODO: implement
  // you have access to all member in cc4s, e.g. cc4s::no, ...
}

