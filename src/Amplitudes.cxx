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
  ai->set_name("Tai");
  {
    int syms[] = {NS, NS, NS, NS};
    abij = new Tensor<>(4,V->get(ABIJ).lens,syms, *V->get(ABIJ).wrld, "Tabij");
    get(ABIJ)["abij"] = V->get(ABIJ)["abij"];
  }
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

