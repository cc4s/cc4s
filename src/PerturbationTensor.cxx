/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <PerturbationTensor.hpp>
#include <Exception.hpp>
#include <iostream>

using namespace cc4s;
using namespace CTF;

/**
 * a,b,c,d,e,f: particle indices
 * g,h: plane wave indicies
 * i,j,k,l,m,n: hole indices
 * x,y: slices of particle indices (in tensors involving more than 2 particles)
 */
Idx_Tensor PerturbationTensor::operator [](char const *indexMap) {
  int mapLength, no, nv, nx, ng;
  mapLength = strlen(indexMap);
  char stdIndexMap[mapLength+1];
  stdIndexMap[mapLength]='\0';
  no = 0;
  nv = 0;
  nx = 0;
  ng = 0;
  for (int i(0); i < mapLength; ++i) {
    if (indexMap[i] > 'z') {
      std::stringstream stream("");
      stream << "Invalid index: " << indexMap;
      throw new Exception(stream.str());
    } else if (indexMap[i] >= 'x') {
      stdIndexMap[i] = 'x' + nx;
      nx++;
    } else if (indexMap[i] >= 'i') {
      stdIndexMap[i] = 'i' + no;
      no++;
    } else if (indexMap[i] >= 'g') {
      stdIndexMap[i] = 'g' + ng;
      ng++;
    } else if (indexMap[i] >= 'a') {
      stdIndexMap[i] = 'a' + nv;
      nv++;
    } else {
      std::stringstream stream("");
      stream << "Invalid index: " << indexMap;
      throw new Exception(stream.str());
    }
  }
  return get(stdIndexMap, indexMap);
}

