/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "Amplitudes.hpp"
#include "Exception.hpp"
#include <iostream>

using namespace CTF;

/**
 * \brief Allocate all tensors
 */
Amplitudes::Amplitudes(
  CoulombIntegrals *V, World *world, Options const * options
): PerturbationTensor(world, options) {
  if (world->rank == 0) std::cout << "Initializing amplitudes...";
  ai = new Tensor<>(V->ai);
  ai->set_name("Tai");
  {
    abij = new Tensor<>(V->abij);
//    int syms[] = {NS, NS, NS, NS};
//    abij = new Tensor<>(4, V->abij->lens, syms, *world, "Tabij");
//    (*abij)["abij"] = (*V)["abij"];
  }
  if (world->rank == 0) std::cout << " OK" << std::endl;
}

Amplitudes::~Amplitudes() {
  delete ai;
  delete abij;
}

Idx_Tensor Amplitudes::get(char const *stdIndexMap, char const *indexMap) {
  if (0 == strcmp(stdIndexMap, "ai")) return (*ai)[indexMap];
  if (0 == strcmp(stdIndexMap, "abij")) return (*abij)[indexMap];
  {
    std::stringstream stream("");
    stream << "Cannot fetch Amplitudes tensor part " << indexMap;
    throw new Exception(stream.str());
  }
}

