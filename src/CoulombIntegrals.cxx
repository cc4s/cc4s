/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <CoulombIntegrals.hpp>
#include <Cc4s.hpp>
#include <util/Log.hpp>
#include <Exception.hpp>
#include <iostream>

using namespace cc4s;
using namespace CTF;

/**
 * \brief Allocate all Coulomb integral tensors and calculate them
 * from the given chi tensors.
 */
CoulombIntegrals::CoulombIntegrals(
  Chi *chiReal, Chi *chiImag
): PerturbationTensor(), chiR(chiReal), chiI(chiImag), v() {
  int nv = chiR->nv;
  int no = chiR->no;
  bool profile = Cc4s::options->profile;
/*
  int symsASAS[] = { AS, NS, AS, NS };
  int symsASNS[] = { AS, NS, NS, NS };
  int symsNSNS[] = { NS, NS, NS, NS };
  int symsNSAS[] = { NS, NS, AS, NS };
*/
  int symsNSNS[] = { NS, NS, NS, NS };
  int vvvv[] = { nv, nv, nv, nv };
  int vvvo[] = { nv, nv, nv, no };
  int vovv[] = { nv, no, nv, nv };
  int vovo[] = { nv, no, nv, no };
  int vvoo[] = { nv, nv, no, no };
  int oovv[] = { no, no, nv, nv };
  int vooo[] = { nv, no, no, no };
  int voov[] = { nv, no, no, nv };
  int oovo[] = { no, no, nv, no };
  int oooo[] = { no, no, no, no };

  // allocate the tensors, assign them to the respective variable
  // and add them to the tensor map for further manipulation
  add(i = new Vector<>(no, *Cc4s::world, "Ei", profile));
  add(a = new Vector<>(nv, *Cc4s::world, "Ea", profile));
  add(ij = new Matrix<>(no, no, NS, *Cc4s::world, "Fij", profile));
  add(ia = new Matrix<>(no, nv, NS, *Cc4s::world, "Fia", profile));
  add(ai = new Matrix<>(nv, no, NS, *Cc4s::world, "Fai", profile));
  add(ab = new Matrix<>(nv, nv, NS, *Cc4s::world, "Fab", profile));
  add(ijkl = new Tensor<>(4, oooo, symsNSNS, *Cc4s::world, "Vijkl", profile));
  add(ijak = new Tensor<>(4, oovo, symsNSNS, *Cc4s::world, "Vijak", profile));
  add(aijk = new Tensor<>(4, vooo, symsNSNS, *Cc4s::world, "Vaijk", profile));
  add(aijb = new Tensor<>(4, voov, symsNSNS, *Cc4s::world, "Vaijb", profile));
  add(ijab = new Tensor<>(4, oovv, symsNSNS, *Cc4s::world, "Vijab", profile));
  add(abij = new Tensor<>(4, vvoo, symsNSNS, *Cc4s::world, "Vabij", profile));
  add(aibj = new Tensor<>(4, vovo, symsNSNS, *Cc4s::world, "Vaibj", profile));
  if (Cc4s::options->storeV) {
    add(aibc = new Tensor<>(4, vovv, symsNSNS, *Cc4s::world, "Vaibc", profile));
    add(abci = new Tensor<>(4, vvvo, symsNSNS, *Cc4s::world, "Vabci", profile));
    add(abcd = new Tensor<>(4, vvvv, symsNSNS, *Cc4s::world, "Vabcd", profile));
  } else {
    aibc = nullptr;
    abci = nullptr;
    abcd = nullptr;
  }
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
  if (t != nullptr) {
    return (*t)[indexMap];
  } else {
    std::stringstream stream("");
    stream << "Cannot fetch CoulombIntegrals tensor part " << indexMap;
    throw new Exception(stream.str());
  }
}

Tensor<> CoulombIntegrals::getSlice(int a, int b) {
  int nv = chiR->nv;
  int no = chiR->no;
  // NOTE: width of sliced hardcoded
  int w = Cc4s::options->nw;
  Tensor<> Rgxc(chiR->getSlice(no+a,no+a+w, no,no+nv)); Rgxc.set_name("Rgxc");
  Tensor<> Rgcy(chiR->getSlice(no,no+nv, no+b,no+b+w)); Rgcy.set_name("Rgcy");
  Tensor<> Igxc(chiI->getSlice(no+a,no+a+w, no,no+nv)); Igxc.set_name("Igxc");
  Tensor<> Igcy(chiI->getSlice(no,no+nv, no+b,no+b+w)); Igcy.set_name("Igcy");
  int lens[] = {Rgxc.lens[1], Rgcy.lens[2], Rgxc.lens[2], Rgcy.lens[1]};
  int syms[] = {NS, NS, NS, NS};
  Tensor<> Vxycd(4, lens, syms, *Cc4s::world, "Vxycd", Cc4s::options->profile);
  Vxycd["xycd"] =  Rgxc["gxc"]*Rgcy["gdy"];
//  Vxycd["xycd"] -= Rgxc["gxd"]*Rgcy["gcy"];
  Vxycd["xycd"] += Igxc["gxc"]*Igcy["gdy"];
//  Vxycd["xycd"] -= Igxc["gxd"]*Igcy["gcy"];
  // NOTE: ctf double counts if lhs tensor is AS
//  Vxycd["xycd"] = 0.5 * Vxycd["xycd"];
  return Vxycd;
}

/* Fetch direct integrals only */
void CoulombIntegrals::fetch(Tensor<> &t, char const *indexMap) {
  if (strlen(indexMap) == 4) {
    LOG(0) << "Calculating V" << indexMap << "...";
    // 4 point tensors:
    char dirL[4] = {'g', indexMap[0], indexMap[2], 0};
    char dirR[4] = {'g', indexMap[3], indexMap[1], 0};
    t[indexMap]  = (*chiR)[dirL] * (*chiR)[dirR];
//    t[indexMap] -= (*chiR)[excL] * (*chiR)[excR];
    t[indexMap] += (*chiI)[dirL] * (*chiI)[dirR];
//    t[indexMap] -= (*chiI)[excL] * (*chiI)[excR];

    int *syms = t[indexMap].parent->sym;
    if (syms[0] == AS || syms[2] == AS) {
      // NOTE: ctf double counts if lhs tensor is AS
      t[indexMap] = 0.5 * t[indexMap];
    }
    LOG(0) << " OK" << std::endl;
  }
}

/**
 * \brief Fetches all tensors elements
 */
void CoulombIntegrals::fetch() {
  // FIXME: how to calculate Vij, Via, Vai, Vab
/*
  chiR->readRandom(ij, 7);
  chiR->readRandom(ia, 8);
  chiR->readRandom(ai, 9);
  chiR->readRandom(ab, 10);
*/

  // TODO: only set fully anti-symmetrized tensor, otherwise
  // we get wrong results

  // calculate the other Coulomb integrals from the chis
  for (auto p(v.begin()); p != v.end(); ++p) {
    fetch(*p->second, p->first.c_str());
  }
}

