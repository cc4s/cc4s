/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "Amplitudes.hpp"
#include <iostream>

using namespace CTF;

/**
 * \brief Allocate all tensors
 */
Amplitudes::Amplitudes(CoulombIntegrals *V) {
  std::cout << "Initializing amplitudes...";
  ai = new Tensor<>(V->get(AI));
  abij = new Tensor<>(V->get(ABIJ));
  std::cout << " OK." << std::endl;
}

Tensor<> &Amplitudes::get(Part part) {
  switch (part) {
    case AI: return *ai;
    case ABIJ: return *abij;
    default: {
      std::stringstream stream("Cannot fetch tensor T part #");
      stream << part;
      throw new Exception(stream.str());
    }
  }
}

