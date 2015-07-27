/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "CoulombIntegrals.hpp"
#include "Exception.hpp"
#include <iostream>

using namespace CTF;

/**
 * \brief Allocate all tensors
 */
CoulombIntegrals::CoulombIntegrals(
  Chi *chiReal_, Chi *chiImag_
): chiReal(chiReal_), chiImag(chiImag_) {
  int nv = chiReal->ai->lens[1];
  int no = chiReal->ai->lens[2];
  World *world = chiReal->ai->wrld;
  bool profile = chiReal->ai->profile;
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

CoulombIntegrals::~CoulombIntegrals() {
  delete abcd;
  delete abij;
  delete ai;
}

Idx_Tensor CoulombIntegrals::get(char const *stdIndexMap, char const *indexMap){
//    printf("indices %s are %s\n",idx_map_,new_idx_map);
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
//  int syms[] = {a == b ? SY : NS, NS, SY, NS};
  int syms[] = {a == b ? NS : NS, NS, NS, NS};
  Tensor<> Vxycd(4, lens, syms, *Rgxc.wrld, "Vxycd", Rgxc.profile);
  Vxycd["xycd"] =  Rgxc["gxc"]*Rgyc["gyd"];
  Vxycd["xycd"] -= Rgxc["gxd"]*Rgyc["gyc"];
  Vxycd["xycd"] += Igxc["gxc"]*Igyc["gyd"];
  Vxycd["xycd"] -= Igxc["gxd"]*Igyc["gyc"];
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

  // calculate the other Coulomb integrals from the chis
  (*abij)["abij"] =  (*chiReal)["gai"]*(*chiReal)["gbj"];
  (*abij)["abij"] -= (*chiReal)["gaj"]*(*chiReal)["gbi"];
  (*abij)["abij"] += (*chiImag)["gai"]*(*chiImag)["gbj"];
  (*abij)["abij"] -= (*chiImag)["gaj"]*(*chiImag)["gbi"];
  // NOTE: only calculate for testing
  (*abcd)["abcd"] =  (*chiReal)["gac"]*(*chiReal)["gbd"];
  (*abcd)["abcd"] -= (*chiReal)["gad"]*(*chiReal)["gbc"];
  (*abcd)["abcd"] += (*chiImag)["gac"]*(*chiImag)["gbd"];
  (*abcd)["abcd"] -= (*chiImag)["gad"]*(*chiImag)["gbc"];
}

