/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Chi.hpp>
#include <Exception.hpp>
#include <Cc4s.hpp>
#include <iostream>
#include <string>
#include <iostream>
#include <fstream>

using namespace cc4s;
using namespace CTF;

Chi::Chi(
  int nG_, int no_, int nv_
): PerturbationTensor(), nG(nG_), no(no_), nv(nv_) {
  int lens[] = {nG, no+nv, no+nv};
  int syms[] = {NS, NS, NS};
  gpq = new Tensor<>(3, lens, syms, *Cc4s::world,"Xgpq",Cc4s::options->profile);
}

Idx_Tensor Chi::get(char const *stdIndexMap, char const *indexMap) {
  // NOTE: PerturbationTensor does not know about p,q,.. indices in stdIndexMap
  if (0 == strcmp(indexMap, "gpq")) return (*gpq)[indexMap];
  if (indexMap[0] == 'g') {
    int start[] = {0, 0, 0};
    int end[] = {gpq->lens[0], 0, 0};
    for (int index = 1; index < 3; ++index) {
      if (stdIndexMap[index] <= 'b') {
        start[index] = no, end[index] = no+nv;
      } else if (stdIndexMap[index] <= 'j') {
        start[index] = 0; end[index] = no;
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

Chiai::Chiai(
  int nG_, int no_, int nv_
): PerturbationTensor(), nG(nG_), no(no_), nv(nv_) {
  int lens[] = {nG, nv, no};
  int syms[] = {NS, NS, NS};
  gai = new Tensor<>(3, lens, syms, *Cc4s::world,"Xgai",Cc4s::options->profile);
}

Idx_Tensor Chiai::get(char const *stdIndexMap, char const *indexMap) {
  // NOTE: PerturbationTensor does not know about p,q,.. indices in stdIndexMap
  if (0 == strcmp(indexMap, "gai")) return (*gai)[indexMap];
  if (indexMap[0] == 'g') {
    int start[] = {0, 0, 0};
    int end[] = {gai->lens[0], 0, 0};
    for (int index = 1; index < 3; ++index) {
      if (stdIndexMap[index] <= 'b') {
        start[index] = no, end[index] = no+nv;
      } else if (stdIndexMap[index] <= 'j') {
        start[index] = 0; end[index] = no;
      }
    }
    //FIXME: memory leak from dynamic allocation of Tensor object
    // find a better solution...
    return (*new Tensor<>(gai->slice(start, end)))[indexMap];
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

Chiai::~Chiai() {
  delete gai;
}
