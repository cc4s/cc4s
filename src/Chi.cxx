/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "Chi.hpp"
#include "Exception.hpp"
#include <iostream>
#include <string>

using namespace CTF;

Chi::Chi(
  World *world, Options const *options, int seed
): PerturbationTensor(world, options) {
  // keep the chi tensors in the memory for now
  {
    int lens[] = {options->nG, options->nv, options->nv};
    int syms[] = {NS, NS, NS};
    gab = new Tensor<>(
      3, lens, syms, *world, "Xgab", options->profile
    );
  }
  {
    int lens[] = {options->nG, options->nv, options->no};
    int syms[] = {NS, NS, NS};
    gai = new Tensor<>(
      3, lens, syms, *world, "Xgai", options->profile
    );
  }
  {
    int lens[] = {options->nG, options->no, options->no};
    int smys[] = {NS, NS, NS};
    gij = new Tensor<>(
      3, lens, smys, *world, "Xgij", options->profile
    );
  }
  readRandom(gab, 0+seed);
  readRandom(gai, 2+seed);
  readRandom(gij, 4+seed);
  // FIXME: symmetrize: symmetries should be treated at reading
  (*gab)["gab"] += (*gab)["gba"];
  (*gij)["gij"] += (*gij)["gji"];
}

Idx_Tensor Chi::get(char const *stdIndexMap, char const *indexMap) {
  if (0 == strcmp(stdIndexMap, "gij")) return (*gij)[indexMap];
  if (0 == strcmp(stdIndexMap, "gai")) return (*gai)[indexMap];
  if (0 == strcmp(stdIndexMap, "gab")) return (*gab)[indexMap];	
  {
    std::stringstream stream("");
    stream << "Cannot fetch Chi tensor part " << stdIndexMap <<
      " with index names " << indexMap;
    throw new Exception(stream.str());
  }
}

Tensor<> Chi::getSlice(int a) {
  int na = std::min(gab->lens[1]-a, gai->lens[2]);
  int start[] = {0, a, 0};
  int end[] = {gab->lens[0], a+na, gab->lens[2]};
  return gab->slice(start, end);
}

Chi::~Chi() {
  delete gab; delete gai; delete gij;
}

void Chi::readRandom(Tensor<> *tensor, int seed) {
  int64_t indicesCount;
  int64_t *indices;
  double *values;
  if (world->rank == 0) {
    std::cout << "Fetching randomly " << tensor->get_name() << " ..." <<
      std::endl;
  }
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

