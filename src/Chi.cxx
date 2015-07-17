/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "Chi.hpp"
#include "Exception.hpp"
#include <iostream>

using namespace CTF;

Chi::Chi(World *world, Options const &options) {
  // keep the chi tensors in the memory for now
  {
    int lens[] = {options.nG, options.nv, options.nv};
    int syms[] = {NS, SY, NS};
    ab = new Tensor<>(3, lens, syms, *world, "Xab", options.profile);
  }
  {
    int lens[] = {options.nG, options.nv, options.no};
    int syms[] = {NS, NS, NS};
    ai = new Tensor<>(3, lens, syms, *world, "Xai", options.profile);
  }
  {
    int lens[] = {options.nG, options.no, options.no};
    int smys[] = {NS, NS, NS};
    ij = new Tensor<>(3, lens, smys, *world, "Xij", options.profile);
  }
  readRandom(ab);
  readRandom(ai);
  readRandom(ij);
}

Tensor<> &Chi::get(ChiPart part) {
  switch (part) {
    case GAB: return *ab;
    case GAI: return *ai;
    case GIJ: return *ij;
    default:
      throw new Exception("Cannot fetch chi tensor part #" + part);
  }
}

Chi::~Chi() {
  delete ab; delete ai; delete ij;
}

void Chi::readRandom(Tensor<> *tensor) {
  int64_t indicesCount;
  int64_t *indices;
  double *values;
  std::cout << "Fetching " << tensor->get_name() << " ..." << std::endl;
  tensor->read_local(&indicesCount, &indices, &values);
  for (int64_t j(0); j < indicesCount; ++j) {
    values[j] = ((indices[j]*13)%13077)/13077. -.5;
  }
  tensor->write(indicesCount, indices, values);
  free(indices); free(values);
}

/**
 * \brief Reads the chi amplitudes from disk
 */
void Chi::read() {
  // TODO: implement
}

