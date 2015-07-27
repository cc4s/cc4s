/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "Chi.hpp"
#include "Exception.hpp"
#include <iostream>
#include <string>

using namespace CTF;

Chi::Chi(World *world, Options const &options, int seed) {
  // keep the chi tensors in the memory for now
  {
    int lens[] = {options.nG, options.nv, options.nv};
    int syms[] = {NS, NS, NS};
    ab = new Tensor<>(
      3, lens, syms, *world, "Xab", options.profile
    );
  }
  {
    int lens[] = {options.nG, options.nv, options.no};
    int syms[] = {NS, NS, NS};
    ai = new Tensor<>(
      3, lens, syms, *world, "Xai", options.profile
    );
  }
  {
    int lens[] = {options.nG, options.no, options.no};
    int smys[] = {NS, SY, NS};
    ij = new Tensor<>(
      3, lens, smys, *world, "Xij", options.profile
    );
  }
  readRandom(ab, 0+seed);
  readRandom(ai, 2+seed);
  readRandom(ij, 4+seed);
}

Idx_Tensor Chi::get(char const *stdIndexMap, char const *indexMap) {
  std::cout << stdIndexMap << ":" << indexMap << std::endl;
  if (0 == strcmp(stdIndexMap, "gij")) return (*ij)[indexMap];
  if (0 == strcmp(stdIndexMap, "gai")) return (*ai)[indexMap];
  if (0 == strcmp(stdIndexMap, "gab")) return (*ab)[indexMap];	
  {
    std::stringstream stream("");
    stream << "Cannot fetch Chi tensor part " << stdIndexMap <<
      " with index names " << indexMap;
    throw new Exception(stream.str());
  }
}

Tensor<> Chi::getSlice(int a) {
  int na = std::min(ab->lens[1]-a, ai->lens[2]);
  int start[] = {0, a, 0};
  int end[] = {ab->lens[0], a+na, ab->lens[2]};
  return ab->slice(start, end);
}

Chi::~Chi() {
  delete ab; delete ai; delete ij;
}

void Chi::readRandom(Tensor<> *tensor, int seed) {
  int64_t indicesCount;
  int64_t *indices;
  double *values;
  std::cout << "Fetching " << tensor->get_name() << " ..." << std::endl;
  tensor->read_local(&indicesCount, &indices, &values);
  for (int64_t j(0); j < indicesCount; ++j) {
    values[j] = ((indices[j]*13 + seed)%13077)/13077. -.5;
    values[j] = ((indices[j]*13 + seed+1)%13077)/13077. -.5;
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

