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
  
  Matrix<> X((int)rank, np, NS, *chiR->wrld, "XRp", chiR->profile);
  Matrix<> gamR((int)rank, nG, NS, *chiR->wrld, "gamRRG", chiR->profile);
  Matrix<> gamI((int)rank, nG, NS, *chiR->wrld, "gamIRG", chiR->profile);

  Tensor<> chi0R(3, chiR->lens, chiR->sym, *chiR->wrld, "chi0RGqr", chiR->profile);
  Tensor<> chi0I(3, chiR->lens, chiR->sym, *chiR->wrld, "chi0IGqr", chiR->profile);
  Tensor<> RR(3, chiR->lens, chiR->sym, *chiR->wrld, "RRGqr", chiR->profile);
  Tensor<> RI(3, chiR->lens, chiR->sym, *chiR->wrld, "RIGqr", chiR->profile);
  do {
    chi0R["Gqr"] = X["Rq"]*X["Rr"]*gamR["RG"];
    chi0I["Gqr"] = X["Rq"]*X["Rr"]*gamI["RG"];
    RR["Gqr"] = chi0R["Gqr"] - (*chiR)["Gqr"];
    RI["Gqr"] = chi0I["Gqr"] - (*chiI)["Gqr"];
  } while (false);
}
