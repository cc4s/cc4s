/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "RalsFtodRankDecomposition.hpp"
#include "Exception.hpp"
#include "util/Log.hpp"
#include "util/ComplexTensor.hpp"
#include "util/RandomTensor.hpp"
#include "util/MathFunctions.hpp"
#include "util/IterativePseudoInverter.hpp"
#include <iostream>
#include <fstream>
#include <random>
#include <limits>

using namespace cc4s;
using namespace CTF;

RalsFtodRankDecomposition::RalsFtodRankDecomposition(
  std::vector<Argument const *> const &argumentList
): Algorithm(argumentList) {
}

RalsFtodRankDecomposition::~RalsFtodRankDecomposition() {
}

void RalsFtodRankDecomposition::run() {
  rank = getIntegerArgument("rank");
  Tensor<double> const *chiR(getTensorArgument("chiR"));
  Tensor<double> const *chiI(getTensorArgument("chiI"));
  double epsilon = getRealArgument("epsilon");
  int nG(chiR->lens[0]);
  int np(chiR->lens[1]);
  LOG(3) << "rank=" << rank << std::endl;
  LOG(3) << "nG=" << nG << ", np=" << np << std::endl;
  chi = new Tensor<complex>(
    3, chiR->lens, chiR->sym, *chiR->wrld, "chiGqr", chiR->profile
  );
  toComplexTensor(*chiR, *chiI, *chi);
  // allocate factor tensors
  x = new Matrix<complex>(
    int(rank), np, NS, *chiR->wrld, "xRp", chiR->profile
  );
  gamma = new Matrix<complex>(
    int(rank), nG, NS, *chiR->wrld, "gamRG", chiR->profile
  );
  setRandomTensor(*x);
  double norm(frobeniusNorm(*x));
  (*x)["Rq"] *= 1.0/norm;
  setRandomTensor(*gamma);
  norm = frobeniusNorm(*gamma);
  (*gamma)["RG"] *= 1.0/norm;
  while (true) {
    double lambda;
    std::ifstream lambdaFile;
    lambdaFile.open("lambda");
    lambdaFile >> lambda;
    lambdaFile.close();
    LOG(2) << "lambda=" << lambda << std::endl;
    fit(lambda);
  }
}

void RalsFtodRankDecomposition::fit(double lambda) {
  fitRals("Gqr", *x,'q', *x,'r', *gamma,'G', lambda);
  fitRals("Gqr", *x,'r', *gamma,'G', *x,'q', lambda);
  fitRals("Gqr", *gamma,'G', *x,'q', *x,'r', lambda);

  int bcLens[] = { x->lens[0], x->lens[1], x->lens[1] };
  int bcSyms[] = { NS, NS, NS };
  Tensor<complex> bc(3, bcLens, bcSyms, *chi->wrld, "bcRjk", chi->profile);
  bc["Rqr"] = (*x)["Rq"] * (*x)["Rr"];
  Tensor<complex> chi0(3, chi->lens, chi->sym, *chi->wrld, "chi0Gqr", chi->profile);
  chi0["Gqr"] = (*gamma)["RG"] * bc["Rqr"];
  chi0["Gqr"] -= (*chi)["Gqr"];
  double R(frobeniusNorm(chi0));
  LOG(2) << R*R << std::endl;
}

void RalsFtodRankDecomposition::fitAls(
  char const *indicesChi,
  Tensor<complex> &b, char const idxB, Tensor<complex> &c, char const idxC,
  Tensor<complex> &a, char const idxA
) {
  Matrix<complex> bb(rank, rank, NS, *chi->wrld, "bbRS", chi->profile);
  Matrix<complex> gramian(rank, rank, NS, *chi->wrld, "gRS", chi->profile);
  Bivar_Function<complex> fDot(&dot<complex>);
  LOG(4) << "building Gramian..." << std::endl;
  bb.contract(1.0,b,"Rj", b,"Sj", 0.0,"SR", fDot);
  gramian.contract(1.0,c,"Rk", c,"Sk", 0.0,"SR", fDot);
  gramian["SR"] *= bb["SR"];
  int bcLens[] = { b.lens[0], b.lens[1], c.lens[1] };
  int bcSyms[] = { NS, NS, NS };
  Tensor<complex> bc(3, bcLens, bcSyms, *chi->wrld, "bcRjk", chi->profile);
  bc["Sjk"] = b["Sj"] * c["Sk"];
  LOG(4) << "inverting Gramian..." << std::endl;
  IterativePseudoInverter<complex> gramianInverter(gramian);
  bc["Rjk"] = bc["Sjk"] * gramianInverter.invert()["SR"];
  char const indicesA[] = { 'R', idxA, 0 };
  char const indicesBC[] = { 'R', idxB, idxC, 0 };
  LOG(4) << "building new estimate: A[" << indicesA <<
    "] = chi[" << indicesChi << "] * bc[" << indicesBC << "] ..." << std::endl;
  a.contract(1.0, *chi,indicesChi, bc,indicesBC, 0.0, indicesA, fDot);
  LOG(4) << "done" << std::endl;
}

void RalsFtodRankDecomposition::fitRals(
  char const *indicesChi,
  Tensor<complex> &b, char const idxB, Tensor<complex> &c, char const idxC,
  Tensor<complex> &a, char const idxA,
  double lambda
) {
  Matrix<complex> bb(rank, rank, NS, *chi->wrld, "bbRS", chi->profile);
  Matrix<complex> gramian(rank, rank, NS, *chi->wrld, "gRS", chi->profile);
  Bivar_Function<complex> fDot(&dot<complex>);
  LOG(4) << "building Gramian..." << std::endl;
  bb.contract(1.0,b,"Rj", b,"Sj", 0.0,"SR", fDot);
  gramian.contract(1.0,c,"Rk", c,"Sk", 0.0,"SR", fDot);
  gramian["SR"] *= bb["SR"];
  gramian["RR"] += lambda;
  int bcLens[] = { b.lens[0], b.lens[1], c.lens[1] };
  int bcSyms[] = { NS, NS, NS };
  Tensor<complex> bc(3, bcLens, bcSyms, *chi->wrld, "bcRjk", chi->profile);
  bc["Sjk"] = b["Sj"] * c["Sk"];
  char const indicesA[] = { 'R', idxA, 0 };
  char const indicesBC[] = { 'R', idxB, idxC, 0 };
  a.contract(1.0, *chi,indicesChi, bc,indicesBC, lambda, indicesA, fDot);

  LOG(4) << "inverting Gramian..." << std::endl;
  IterativePseudoInverter<complex> gramianInverter(gramian);
  a.contract(1.0, a,"Si", gramianInverter.invert(), "SR", 0.0, "Ri", fDot);
}

void RalsFtodRankDecomposition::test(CTF::World *world) {
  Matrix<complex> t(5,7, NS, *world, "Tij");
  Matrix<complex> a(5,10, NS, *world, "AiR");
  Matrix<complex> b(7,10, NS, *world, "BjR");
  setRandomTensor(t);
  setRandomTensor(b);
  Bivar_Function<complex> fDot(&dot<complex>);
  Matrix<complex> gramian(10,10, NS, *world, "GSR");
  gramian.contract(1.0,b,"jR", b,"jS", 0.0,"SR", fDot);
  Matrix<complex> bg(7,10, NS, *world, "BjR");
  IterativePseudoInverter<complex> gramianInverter(gramian);
  bg["jR"] = b["jS"] * gramianInverter.invert()["SR"];
  a.contract(1.0, t,"ij", bg,"jR", 0.0,"iR", fDot);

  Matrix<complex> ab(5,7, NS, *world, "Tij");
  ab["ij"] = a["iR"] * b["jR"];
  ab["ij"] -= t["ij"];
  logMatrix(2, ab);
  double error(frobeniusNorm(ab));
  LOG(2) << "error=" << error << std::endl;
}

