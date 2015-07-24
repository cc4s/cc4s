/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "CoulombIntegrals.hpp"
#include "Exception.hpp"
#include <iostream>

using namespace CTF;

/**
 * \brief Allocate all tensors
 */
CoulombIntegrals::CoulombIntegrals(Chi *chi_): chi(chi_) {
  int nv = chi->get(GAI).lens[1];
  int no = chi->get(GAI).lens[2];
  World *world = chi->get(GAI).wrld;
  bool profile = chi->get(GAI).profile;
  a = new Vector<>(nv, *world, "Va", profile);
  i = new Vector<>(no, *world, "Vi", profile);
  {
    int lens[] = {nv, no};
    int syms[] = {NS, NS};
    ai = new Tensor<>(2, lens, syms, *world, "Vai", profile);
  }
  {
    int lens[] = {nv, nv, no, no};
    int syms[] = {NS, NS, NS, NS};
    abij = new Tensor<>(4, lens, syms, *world, "Vabij",profile);
  }
// NOTE: only for testing
  {
    int lens[] = {nv, nv, nv, nv};
//    int syms[] = {SY, NS, SY, NS};
    int syms[] = {NS, NS, NS, NS};
    abcd = new Tensor<>(4, lens, syms, *world, "Vabcd",profile);
  }
  fetch();
}

Tensor<> &CoulombIntegrals::get(Part part) {
  switch (part) {
    case A: return *a;
    case I: return *i;
    case AI: return *ai;
    case ABIJ: return *abij;
    case ABCD:
      return *abcd;
      throw new Exception("Cannot fetch entire ABCD tensor in memory");
    default: {
      std::stringstream stream("Cannot fetch tensor V part #");
      stream << part;
      throw new Exception(stream.str());
    }
  }
}

Tensor<> CoulombIntegrals::getSlice(Part part, int a, int b) {
  switch (part) {
    case ABCD: {
      Tensor<> Xgxc(chi->getSlice(GAB, a));
      Tensor<> Xgyd(chi->getSlice(GAB, b));
      int lens[] = {Xgxc.lens[1], Xgyd.lens[1], Xgxc.lens[2], Xgyd.lens[2]};
//      int syms[] = {a == b ? SY : NS, NS, SY, NS};
      int syms[] = {a == b ? NS : NS, NS, NS, NS};
      Tensor<> Vxycd(4, lens, syms, *Xgxc.wrld, "Vxycd", Xgxc.profile);
      Vxycd["xycd"] = Xgxc["gxc"] * Xgyd["gyd"];
      return Vxycd;
    }
    default: {
      std::stringstream stream("Cannot fetch slice of tensor V part #");
      stream << part;
      throw new Exception(stream.str());
    }
  }
}

/**
 * \brief Fetches all tensors elements
 */
void CoulombIntegrals::fetch() {
  // TODO: read epsilons
  chi->readRandom(a, 3);
  chi->readRandom(i, 4);

  // FIXME: how to calculate Via
  chi->readRandom(ai, 5);

  // calculate the other Coulomb integrals from the chis
  // TODO: use complex numbers and conjugate
  get(ABIJ)["abij"] =
    chi->get(GAI)["gai"] * chi->get(GAI)["gbj"] -
    chi->get(GAI)["gaj"] * chi->get(GAI)["gbi"];
  // NOTE: only calculate for testing
  get(ABCD)["abcd"] = chi->get(GAB)["gac"] * chi->get(GAB)["gbd"];
}

/**
 * \deprecated
 * \brief Calculates a slice xycd of the full tensor abcd
 */
void CoulombIntegrals::calculate_xycd(Tensor<> &xycd, int a, int b) {
  int64_t indicesCount, *indices;
  double *values;
  int nv = xycd.lens[2];
  // lengths in the full Vabcd tensor
  int absoluteLengths[] = {nv, nv, nv, nv};
  int positions[4];
  int64_t absoluteIndex;
  xycd.read_local(&indicesCount, &indices, &values);
  for (int j(0); j < indicesCount; ++j) {
    // get position within slice 
    from_index(xycd, indices[j], positions);
    // compute absolute position in Vabcd
    positions[0] += a;
    positions[1] += b;
    // get index of absolute position
    absoluteIndex = to_index(4, absoluteLengths, positions);
    // use it for calculating the respective element
    // NOTE: Vabcd is the array #4:
    values[j] = ((absoluteIndex*16+4)%13077)/13077. -.5;
  }
  xycd.write(indicesCount, indices, values);
  free(indices), free(values);
}

