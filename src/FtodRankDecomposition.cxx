/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "FtodRankDecomposition.hpp"
#include "Exception.hpp"
#include <iostream>
#include <random>

using namespace CTF;

FtodRankDecomposition::FtodRankDecomposition(
  std::vector<Argument const *> const &argumentList
): Algorithm(argumentList) {
}

FtodRankDecomposition::~FtodRankDecomposition() {
}

void FtodRankDecomposition::run() {
  rank = getIntegerArgument("rank");
  chiR = const_cast<Tensor<> *>(getTensorArgument("chiR"));
  chiI = const_cast<Tensor<> *>(getTensorArgument("chiI"));
  std::cout << "rank=" << rank << std::endl;
  int nG(chiR->lens[0]);
  int np(chiR->lens[1]);
  std::cout << "nG=" << nG << ", np=" << np << std::endl;
  // allocate decompisiton tensors
  X = new Matrix<>(int(rank), np, NS, *chiR->wrld, "XRp", chiR->profile);
  gamR = new Matrix<>(int(rank), nG, NS, *chiR->wrld, "gamRRG", chiR->profile);
  gamI = new Matrix<>(int(rank), nG, NS, *chiR->wrld, "gamIRG", chiR->profile);

  // allocate intermediate tensors:
  // tensor rank approximation
  chi0R = new Tensor<>(3, chiR->lens, chiR->sym, *chiR->wrld, "chi0RGqr", chiR->profile);
  chi0I = new Tensor<>(3, chiR->lens, chiR->sym, *chiR->wrld, "chi0IGqr", chiR->profile);
  // residuum
  RR = new Tensor<>(3, chiR->lens, chiR->sym, *chiR->wrld, "RRGqr", chiR->profile);
  RI = new Tensor<>(3, chiR->lens, chiR->sym, *chiR->wrld, "RIGqr", chiR->profile);

  // gradients
  dX = new Matrix<>(int(rank), np, NS, *chiR->wrld, "dXRp", chiR->profile);
  dGamR = new Matrix<>(int(rank), nG, NS, *chiR->wrld, "dGamRRG", chiR->profile);
  dGamI = new Matrix<>(int(rank), nG, NS, *chiR->wrld, "dGamIRG", chiR->profile);

  initializeX();
  initializeGam();
  testGradient();
  sX = new Matrix<>(*dX);
  sGamR = new Matrix<>(*dGamR);
  sGamI = new Matrix<>(*dGamI);
  lineSearchX();
}


void FtodRankDecomposition::calculateChi0() {
  // TODO: try to keep this tensor allocated
  int np(chiR->lens[1]);
  int xLens[] = {int(rank), np, np};
  int xSyms[] = {NS, NS, NS};
  Tensor<> x(3, xLens, xSyms, *chiR->wrld, "xRqr", chiR->profile);

  x["Rqr"] = (*X)["Rq"]*(*X)["Rr"];
  (*chi0R)["Gqr"] = x["Rqr"]*(*gamR)["RG"];
  x["Rqr"] = (*X)["Rq"]*(*X)["Rr"];
  (*chi0I)["Gqr"] = x["Rqr"]*(*gamI)["RG"];
}

// TODO: find proper place for it
double sqr(double x) { return x*x; }

void FtodRankDecomposition::calculateResiduum() {
  // calculate the residuum
  (*RR)["Gqr"] = (*chi0R)["Gqr"] - (*chiR)["Gqr"];
  (*RI)["Gqr"] = (*chi0I)["Gqr"] - (*chiI)["Gqr"];
  R = sqr(RR->norm2()) + sqr(RI->norm2());
}

void FtodRankDecomposition::calculateGradient() {
  // TODO: try to keep these tensors allocated
  int nG(chiR->lens[0]);
  int np(chiR->lens[1]);
  int xLens[] = {int(rank), np, np};
  int xSyms[] = {NS, NS, NS};
  Tensor<> x(3, xLens, xSyms, *chiR->wrld, "xRqr", chiR->profile);
  int gLens[] = {int(rank), nG, np};
  int gSyms[] = {NS, NS, NS};
  Tensor<> g(3, gLens, gSyms, *chiR->wrld, "gRGq", chiR->profile);

  // calculate the gradient of X
  x["Rqr"] = (*RR)["Gqr"]*(*gamR)["RG"];
  // calculate 2*(x["Rpr"]*(*X)["Rr"] + x["Rpq"]*(*X)["Rq"])
  (*dX)["Rp"] = 4.0*x["Rpr"]*(*X)["Rr"];
  x["Rqr"] = (*RI)["Gqr"]*(*gamI)["RG"];
  (*dX)["Rp"] += 4.0*x["Rpr"]*(*X)["Rr"];

  // calculate the gradient of gamR and gamI
  g["RGq"] = (*RR)["Gqr"]*(*X)["Rr"];
  (*dGamR)["RG"] = 2.0*g["RGq"]*(*X)["Rq"];
  g["RGq"] = (*RI)["Gqr"]*(*X)["Rr"];
  (*dGamI)["RG"] = 2.0*g["RGq"]*(*X)["Rq"];
}

void FtodRankDecomposition::initializeRandom(Tensor<> &t, int64_t seed) {
  int64_t indicesCount, *indices;
  double *values;
  std::mt19937 random;
  random.seed(seed);
  std::normal_distribution<double> normalDistribution(0.0, 1.0);
  t.read_local(&indicesCount, &indices, &values);
  for (int64_t i(0); i < indicesCount; ++i) {
    values[i] = normalDistribution(random);
  }
  t.write(indicesCount, indices, values);
  free(indices); free(values);
}

void FtodRankDecomposition::initializeX() {
  initializeRandom(*X, chiR->wrld->rank);
}

void FtodRankDecomposition::initializeGam() {
  initializeRandom(*gamR, 2*chiR->wrld->rank+0);
  initializeRandom(*gamI, 2*chiR->wrld->rank+1);
}

void FtodRankDecomposition::lineSearchXPart(
  Tensor<> &chi0, Tensor<> &chi, Tensor<> &gam
) {
  // TODO: try to keep these tensors allocated
  int np(chi.lens[1]);
  // FIXME: no copy needed for the following two allocations
  Tensor<> dChi0(chi);
  Tensor<> ddChi0(chi);
  Tensor<> deltaChi0(chi0);
  deltaChi0["Gqr"] -= chi["Gqr"];
  int xLens[] = {int(rank), np, np};
  int xSyms[] = {NS, NS, NS};
  Tensor<> x(3, xLens, xSyms, *chi.wrld, "xRqr", chi.profile);
  x["Rqr"] = (*sX)["Rq"] * (*X)["Rr"];
  dChi0["Gqr"] = x["Rqr"] * gam["RG"];
  dChi0["Gqr"] += dChi0["Grq"];
  x["Rqr"] = (*sX)["Rq"] * (*sX)["Rr"];
  ddChi0["Gqr"] = x["Rqr"] * gam["RG"];
  Scalar<> s(*chi.wrld);
  s[""] = 2.0*dChi0["Gqr"] * deltaChi0["Gqr"];
  a1 += s.get_val();
  s[""] = 2.0*ddChi0["Gqr"] * deltaChi0["Gqr"];
  s[""] += dChi0["Gqr"] * dChi0["Gqr"];
  a2 += s.get_val();
  s[""] = 2.0*ddChi0["Gqr"] * dChi0["Gqr"];
  a3 += s.get_val();
  s[""] = ddChi0["Gqr"] * ddChi0["Gqr"];
  a4 += s.get_val();
}

double FtodRankDecomposition::lineSearchX() {
  a1 = a2 = a3 = a4 = 0.0;
  lineSearchXPart(*chi0R, *chiR, *gamR);
  lineSearchXPart(*chi0I, *chiI, *gamI);
  std::cout << "a1=" << a1 << std::endl;
  std::cout << "a2=" << a2 << std::endl;
  std::cout << "a3=" << a3 << std::endl;
  std::cout << "a4=" << a4 << std::endl;
  return 0.0;
}

double FtodRankDecomposition::lineSearchGam() {
  return 0.0;
}

void FtodRankDecomposition::optimizeX() {
}

void FtodRankDecomposition::optimizeGam() {
}

void FtodRankDecomposition::testGradient() {
  calculateChi0();
  calculateResiduum();
  calculateGradient();
  double R0(R);
  Matrix<> oldX(*X), oldGamR(*gamR), oldGamI(*gamI);
  double delta(1e-10);
  double accuracy(1e-6);
//  for (double delta(width/100); delta <= width; delta += width/100) {
    (*X)["Rq"] = oldX["Rq"]+delta*(*dX)["Rq"];
    (*gamR)["RG"] = oldGamR["RG"]+delta*(*dGamR)["RG"];
    (*gamI)["RG"] = oldGamI["RG"]+delta*(*dGamI)["RG"];
    calculateChi0();
    calculateResiduum();
    double epsilon(R - R0);
    double dotEpsilon(
      delta*(
        sqr(dX->norm2()) + sqr(dGamR->norm2()) + sqr(dGamI->norm2())
      )
    );
    if (abs(dotEpsilon-epsilon)/epsilon > accuracy)
      throw new Exception("wrong gradient");
//  }
}
