/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "FtodRankDecomposition.hpp"
#include <iostream>

using namespace CTF;

FtodRankDecomposition::FtodRankDecomposition(
  std::vector<Argument const *> const &argumentList
): Algorithm(argumentList) {
  rank = getIntegerArgument("rank");
  chiR = const_cast<Tensor<> *>(getTensorArgument("chiR"));
  chiI = const_cast<Tensor<> *>(getTensorArgument("chiI"));
  std::cout << "rank=" << rank << std::endl;
}

FtodRankDecomposition::~FtodRankDecomposition() {
}

void FtodRankDecomposition::run() {
  int nG = chiR->lens[0];
  int np = chiR->lens[1];
  std::cout << "nG=" << nG << ", np=" << np << std::endl;
  // allocate decompisiton tensors
  X = new Matrix<>(int(rank), np, NS, *chiR->wrld, "XRp", chiR->profile);
  gamR = new Matrix<>(int(rank), nG, NS, *chiR->wrld, "gamRRG", chiR->profile);
  gamI = new Matrix<>(int(rank), nG, NS, *chiR->wrld, "gamIRG", chiR->profile);


  // allocate intermediate tensors:
  // gradients
  Matrix<> dX(int(rank), np, NS, *chiR->wrld, "dXRp", chiR->profile);
  Matrix<> dGamR(int(rank), nG, NS, *chiR->wrld, "dGamRRG", chiR->profile);
  Matrix<> dGamI(int(rank), nG, NS, *chiR->wrld, "dGamIRG", chiR->profile);

  Tensor<> chi0R(3, chiR->lens, chiR->sym, *chiR->wrld, "chi0RGqr", chiR->profile);
  Tensor<> chi0I(3, chiR->lens, chiR->sym, *chiR->wrld, "chi0IGqr", chiR->profile);
  Tensor<> RR(3, chiR->lens, chiR->sym, *chiR->wrld, "RRGqr", chiR->profile);
  Tensor<> RI(3, chiR->lens, chiR->sym, *chiR->wrld, "RIGqr", chiR->profile);
  int xLens[] = {int(rank), np, np};
  int xSyms[] = {NS, NS, NS};
  Tensor<> x(3, xLens, xSyms, *chiR->wrld, "xRqr", chiR->profile);
  int gLens[] = {int(rank), nG, np};
  int gSyms[] = {NS, NS, NS};
  Tensor<> g(3, gLens, gSyms, *chiR->wrld, "gRGq", chiR->profile);
  do {
    // calculate the tensor rank decompisiton chi0
    x["Rqr"] = (*X)["Rq"]*(*X)["Rr"];
    chi0R["Gqr"] = x["Rqr"]*(*gamR)["RG"];
    x["Rqr"] = (*X)["Rq"]*(*X)["Rr"];
    chi0I["Gqr"] = x["Rqr"]*(*gamI)["RG"];
    // calculate the residuum
    RR["Gqr"] = chi0R["Gqr"] - (*chiR)["Gqr"];
    RI["Gqr"] = chi0I["Gqr"] - (*chiI)["Gqr"];

    // calculate the gradient of X
    x["Rqr"] = RR["Gqr"]*(*gamR)["RG"];
    (*X)["Rq"] = 2.0*x["Rqr"]*(*X)["Rr"];
    x["Rqr"] = RI["Gqr"]*(*gamI)["RG"];
    (*X)["Rq"] += 2.0*x["Rqr"]*(*X)["Rr"];
    // calculate the gradient of gamR and gamI
    g["RGq"] = RR["Gqr"]*(*X)["Rr"];
    dGamR["RG"] = 2.0*g["RGq"]*(*X)["Rq"];
    g["RGq"] = RI["Gqr"]*(*X)["Rr"];
    dGamI["RG"] = 2.0*g["RGq"]*(*X)["Rq"];
  } while (false);
}

void FtodRankDecomposition::calculateGradient() {
}
