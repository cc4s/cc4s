/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "CoulombIntegrals.hpp"
#include "Exception.hpp"
#include <iostream>

using namespace CTF;

/**
 * \brief Allocate all tensors
 */
CoulombIntegrals::CoulombIntegrals(
  Chi *chiReal_, Chi *chiImag_, World *world, Options const *options
): PerturbationTensor(world, options), chiReal(chiReal_), chiImag(chiImag_) {
  int nv = options->nv;
  int no = options->no;
  a = new Vector<>(nv, *world, "Va", options->profile);
  i = new Vector<>(no, *world, "Vi", options->profile);
  {
    int lens[] = {nv, no};
    int syms[] = {NS, NS};
    ai = new Tensor<>(2, lens, syms, *world, "Vai", options->profile);
  }
  {
    int lens[] = {nv, nv, no, no};
    int syms[] = {NS, NS, AS, NS};
    abij = new Tensor<>(4, lens, syms, *world, "Vabij", options->profile);
  }
  if (options->storeV) {
    int lens[] = {nv, nv, nv, nv};
    int syms[] = {NS, NS, AS, NS};
    abcd = new Tensor<>(4, lens, syms, *world, "Vabcd", options->profile);
  } else {
    abcd = NULL;
  }
  fetch();
}

CoulombIntegrals::~CoulombIntegrals() {
  if (abcd) delete abcd;
  delete abij;
  delete ai;
}

Idx_Tensor CoulombIntegrals::get(char const *stdIndexMap, char const *indexMap){
  if (0 == strcmp("a", stdIndexMap)) return (*a)[indexMap];
  if (0 == strcmp("i", stdIndexMap)) return (*i)[indexMap];
//  if (0 == strcmp("ab", stdIndexMap)) return (*ab)[indexMap];
  if (0 == strcmp("ai", stdIndexMap)) return (*ai)[indexMap];
//  if (0 == strcmp("ia", stdIndexMap)) return (*ia)[indexMap];
//  if (0 == strcmp("ij", stdIndexMap)) return (*ij)[indexMap];
  if (0 == strcmp("abcd", stdIndexMap)) return (*abcd)[indexMap];
//  if (0 == strcmp("abci", stdIndexMap)) return (*abci)[indexMap];
//  if (0 == strcmp("aibc", stdIndexMap)) return (*aibc)[indexMap];
//  if (0 == strcmp("aibj", stdIndexMap)) return (*aibj)[indexMap];
  if (0 == strcmp("abij", stdIndexMap)) return (*abij)[indexMap];
//  if (0 == strcmp("ijab", stdIndexMap)) return (*ijab)[indexMap];
//  if (0 == strcmp("aijk", stdIndexMap)) return (*aijk)[indexMap];
//  if (0 == strcmp("ijak", stdIndexMap)) return (*ijak)[indexMap];
//  if (0 == strcmp("ijkl", stdIndexMap)) return (*ijkl)[indexMap];
  {
    std::stringstream stream("");
    stream << "Cannot fetch CoulombIntegrals tensor part " << indexMap;
    throw new Exception(stream.str());
  }
}

Tensor<> CoulombIntegrals::getSlice(int a, int b) {
  Tensor<> Rgxc(chiReal->getSlice(a)); Rgxc.set_name("Rgxc");
  Tensor<> Rgyc(chiReal->getSlice(b)); Rgyc.set_name("Rgyc");
  Tensor<> Igxc(chiImag->getSlice(a)); Igxc.set_name("Igxc");
  Tensor<> Igyc(chiImag->getSlice(b)); Igyc.set_name("Igyc");
  int lens[] = {Rgxc.lens[1], Rgyc.lens[1], Rgxc.lens[2], Rgyc.lens[2]};
  int syms[] = {a == b ? NS : NS, NS, AS, NS};
  Tensor<> Vxycd(4, lens, syms, *world, "Vxycd", options->profile);
  Vxycd["xycd"] =  Rgxc["gxc"]*Rgyc["gyd"];
  Vxycd["xycd"] -= Rgxc["gxd"]*Rgyc["gyc"];
  Vxycd["xycd"] += Igxc["gxc"]*Igyc["gyd"];
  Vxycd["xycd"] -= Igxc["gxd"]*Igyc["gyc"];
  // NOTE: ctf double counts if lhs tensor is AS
  Vxycd["xycd"] = 0.5 * Vxycd["xycd"];
  return Vxycd;
}

/**
 * \brief Fetches all tensors elements
 */
void CoulombIntegrals::fetch() {
  // TODO: read epsilons
  chiReal->readRandom(a, 5);
  chiReal->readRandom(i, 6);

  // FIXME: how to calculate Via
  chiReal->readRandom(ai, 7);

  // TODO: only set fully anti-symmetrized tensor, otherwise
  // we get wrong results
  // calculate the other Coulomb integrals from the chis
  // TODO: only enter fully anti-symmetrized tensors if AS
  (*abij)["abij"] =  (*chiReal)["gai"]*(*chiReal)["gbj"];
  (*abij)["abij"] -= (*chiReal)["gaj"]*(*chiReal)["gbi"];
  (*abij)["abij"] += (*chiImag)["gai"]*(*chiImag)["gbj"];
  (*abij)["abij"] -= (*chiImag)["gaj"]*(*chiImag)["gbi"];
  // NOTE: ctf double counts if lhs tensor is AS
  (*abij)["abij"] = 0.5 * (*abij)["abij"];
  if (abcd) {
    (*abcd)["abcd"] =  (*chiReal)["gac"]*(*chiReal)["gbd"];
    (*abcd)["abcd"] -= (*chiReal)["gad"]*(*chiReal)["gbc"];
    (*abcd)["abcd"] += (*chiImag)["gac"]*(*chiImag)["gbd"];
    (*abcd)["abcd"] -= (*chiImag)["gad"]*(*chiImag)["gbc"];
    (*abcd)["abcd"] = 0.5 * (*abcd)["abcd"];
    // NOTE: ctf double counts if lhs tensor is AS
  }
}

