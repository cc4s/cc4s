/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "CoulombIntegrals.hpp"
#include "Exception.hpp"
#include <iostream>

using namespace CTF;

/**
 * \brief Allocate all Coulomb integral tensors and calculate them
 * from the given chi tensors.
 */
CoulombIntegrals::CoulombIntegrals(
  Chi *chiReal, Chi *chiImag, World *world, Options const *options
): PerturbationTensor(world, options), chiR(chiReal), chiI(chiImag), v() {
  int nv = options->nv;
  int no = options->no;
  int symsASAS[] = { AS, NS, AS, NS };
  int symsASNS[] = { AS, NS, NS, NS };
  int symsNSNS[] = { NS, NS, NS, NS };
  int symsNSAS[] = { NS, NS, AS, NS };
  int vvvv[] = { nv, nv, nv, nv };
  int vvvo[] = { nv, nv, nv, no };
  int vovv[] = { nv, no, nv, nv };
  int vovo[] = { nv, no, nv, no };
  int vvoo[] = { nv, nv, no, no };
  int oovv[] = { no, no, nv, nv };
  int vooo[] = { nv, no, no, no };
  int oovo[] = { no, no, nv, no };
  int oooo[] = { no, no, no, no };

  // allocate the tensors, assign them to the respective variable
  // and add them to the tensor map for further manipulation
  add(i = new Vector<>(no, *world, "Ei", options->profile));
  add(a = new Vector<>(nv, *world, "Ea", options->profile));
  add(ij = new Matrix<>(no, no, AS, *world, "Fij", options->profile));
  add(ia = new Matrix<>(no, nv, NS, *world, "Fia", options->profile));
  add(ai = new Matrix<>(nv, no, NS, *world, "Fai", options->profile));
  add(ab = new Matrix<>(nv, nv, AS, *world, "Fab", options->profile));
  add(ijkl = new Tensor<>(4, oooo, symsASAS, *world, "Vijkl",options->profile));
  add(ijak = new Tensor<>(4, oovo, symsASNS, *world, "Vijak",options->profile));
  add(aijk = new Tensor<>(4, vooo, symsNSAS, *world, "Vaijk",options->profile));
  add(ijab = new Tensor<>(4, oovv, symsASAS, *world, "Vijab",options->profile));
  add(abij = new Tensor<>(4, vvoo, symsASAS, *world, "Vabij",options->profile));
  add(aibj = new Tensor<>(4, vovo, symsNSNS, *world, "Vaibj",options->profile));
  if (options->storeV) {
    add(aibc = new Tensor<>(4, vovv, symsNSAS,*world,"Vaibc",options->profile));
    add(abci = new Tensor<>(4, vvvo, symsASNS,*world,"Vabci",options->profile));
    add(abcd = new Tensor<>(4, vvvv, symsASAS,*world,"Vabcd",options->profile));
  } else {
    aibc = NULL;
    abci = NULL;
    abcd = NULL;
  }

  // fetch all allocated tensors
  fetch();
}

CoulombIntegrals::~CoulombIntegrals() {
  // iterate through all pairs of the map
  for (auto p(v.begin()); p != v.end(); ++p) {
    // the tensor is the "second" member of the pair
    delete p->second;
  }
}

void CoulombIntegrals::add(Tensor<> *t) {
  // the last characters of the name must be the standard index map
  char const *stdIndexMap = &t->name[1];
  // enter the given tensor into the map for later use
  v[stdIndexMap] = t;
}

Idx_Tensor CoulombIntegrals::get(char const *stdIndexMap, char const *indexMap){
  Tensor<> * t(v[stdIndexMap]);
  if (t != NULL) {
    return (*t)[indexMap];
  } else {
    std::stringstream stream("");
    stream << "Cannot fetch CoulombIntegrals tensor part " << indexMap;
    throw new Exception(stream.str());
  }
}

Tensor<> CoulombIntegrals::getSlice(int a, int b) {
  Tensor<> Rgxc(chiR->getSlice(a)); Rgxc.set_name("Rgxc");
  Tensor<> Rgyc(chiR->getSlice(b)); Rgyc.set_name("Rgyc");
  Tensor<> Igxc(chiI->getSlice(a)); Igxc.set_name("Igxc");
  Tensor<> Igyc(chiI->getSlice(b)); Igyc.set_name("Igyc");
  int lens[] = {Rgxc.lens[1], Rgyc.lens[1], Rgxc.lens[2], Rgyc.lens[2]};
  int syms[] = {a == b ? AS : NS, NS, AS, NS};
  Tensor<> Vxycd(4, lens, syms, *world, "Vxycd", options->profile);
  Vxycd["xycd"] =  Rgxc["gxc"]*Rgyc["gyd"];
  Vxycd["xycd"] -= Rgxc["gxd"]*Rgyc["gyc"];
  Vxycd["xycd"] += Igxc["gxc"]*Igyc["gyd"];
  Vxycd["xycd"] -= Igxc["gxd"]*Igyc["gyc"];
  // NOTE: ctf double counts if lhs tensor is AS
  Vxycd["xycd"] = 0.5 * Vxycd["xycd"];
  return Vxycd;
}

void CoulombIntegrals::fetch(Tensor<> *t, char const *indexMap) {
  if (strlen(indexMap) == 4) {
    if (world->rank == 0) std::cout << "Calculating V" << indexMap << "...";
    // 4 point tensors:
    char dirL[4] = {'g', indexMap[0], indexMap[2], 0 };
    char dirR[4] = {'g', indexMap[1], indexMap[3], 0 };
    char excL[4] = {'g', indexMap[0], indexMap[3], 0 };
    char excR[4] = {'g', indexMap[1], indexMap[2], 0 };
    (*this)[indexMap]  = (*chiR)[dirL] * (*chiR)[dirR];
    (*this)[indexMap] -= (*chiR)[excL] * (*chiR)[excR];
    (*this)[indexMap] += (*chiI)[dirL] * (*chiI)[dirR];
    (*this)[indexMap] -= (*chiI)[excL] * (*chiI)[excR];

    int *syms = (*this)[indexMap].parent->sym;
    if (syms[0] == AS || syms[2] == AS) {
      // NOTE: ctf double counts if lhs tensor is AS
      (*this)[indexMap] = 0.5 * (*this)[indexMap];
    }
    if (world->rank == 0) std::cout << " OK" << std::endl;
  }
}

/**
 * \brief Fetches all tensors elements
 */
void CoulombIntegrals::fetch() {
  // TODO: read epsilons
  chiR->readRandom(a, 5);
  chiR->readRandom(i, 6);

  // FIXME: how to calculate Vij, Via, Vai, Vab
  chiR->readRandom(ij, 7);
  chiR->readRandom(ia, 8);
  chiR->readRandom(ai, 9);
  chiR->readRandom(ab, 10);

  // TODO: only set fully anti-symmetrized tensor, otherwise
  // we get wrong results

  // calculate the other Coulomb integrals from the chis
  for (auto p(v.begin()); p != v.end(); ++p) {
    fetch(p->second, p->first.c_str());
  }
}

