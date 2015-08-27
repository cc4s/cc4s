/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "Chi.hpp"
#include "Exception.hpp"
#include <iostream>
#include <string>
#include <iostream>
#include <fstream>

using namespace CTF;

Chi::Chi(
  World *world, Options const *options
): PerturbationTensor(world, options) {
  int lens[] = {options->nG, options->no+options->nv, options->no+options->nv};
  int syms[] = {NS, NS, NS};
  gpq = new Tensor<>(3, lens, syms, *world, "Xgpq", options->profile);
}

Idx_Tensor Chi::get(char const *stdIndexMap, char const *indexMap) {
  // NOTE: PerturbationTensor does not know about p,q,.. indices in stdIndexMap
  if (0 == strcmp(indexMap, "gpq")) return (*gpq)[indexMap];
  if (indexMap[0] == 'g') {
    int start[] = {0, 0, 0};
    int end[] = {gpq->lens[0], 0, 0};
    for (int index = 1; index < 3; ++index) {
      if (stdIndexMap[index] <= 'b') {
        start[index] = options->no, end[index] = options->no+options->nv;
      } else if (stdIndexMap[index] <= 'j') {
        start[index] = 0; end[index] = options->no;
      }
    }
    //FIXME: memory leak from dynamic allocation of Tensor object
    // find a better solution...
    return (*new Tensor<>(gpq->slice(start, end)))[indexMap];
  } else {
    std::stringstream stream("");
    stream << "Cannot fetch Chi tensor part " << stdIndexMap <<
      " with index names " << indexMap;
    throw new Exception(stream.str());
  }
}

Tensor<> Chi::getSlice(int pStart, int pEnd, int qStart, int qEnd) {
  int start[] = {0, pStart, qStart};
  int end[] = {
    gpq->lens[0], std::min(pEnd, gpq->lens[1]), std::min(qEnd, gpq->lens[2])
  };
  return gpq->slice(start, end);
}

Chi::~Chi() {
  delete gpq;
}

void Chi::readRandom(Tensor<> *tensor, int seed) {
  int64_t indicesCount;
  int64_t *indices;
  double *values;
  if (world->rank == 0) {
    std::cout << "Fetching randomly " << tensor->get_name() << "...";
  }
  tensor->read_local(&indicesCount, &indices, &values);
  for (int64_t j(0); j < indicesCount; ++j) {
    values[j] = ((indices[j]*13 + seed)%13077)/13077. -.5;
    values[j] = ((indices[j]*13 + seed+1)%13077)/13077. -.5;
  }
  tensor->write(indicesCount, indices, values);
  free(indices); free(values);
  if (world->rank == 0) std::cout << " OK" << std::endl;
}
