/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <RalsFtodRankDecomposition.hpp>
#include <Exception.hpp>
#include <util/Log.hpp>
#include <util/ComplexTensor.hpp>
#include <util/RandomTensor.hpp>
#include <util/MathFunctions.hpp>
#include <util/IterativePseudoInverse.hpp>
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
//  double epsilon = getRealArgument("epsilon");
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
  chi0 = new Tensor<complex>(3, chi->lens, chi->sym, *chi->wrld, "chi0Gqr", chi->profile);

  while (true) {
    double lambda;
    std::ifstream lambdaFile;
    lambdaFile.open("lambda");
    lambdaFile >> lambda;
    lambdaFile.close();
    if (lambda < 0.0) return;
    LOG(2) << "lambda=" << lambda << std::endl;
    fit(lambda);
  }
}

void RalsFtodRankDecomposition::fit(double lambda) {
  double deltaG(fitRals("Gqr", *x,'q', *x,'r', *gamma,'G', lambda));
  double deltaX1(fitRals("Gqr", *x,'r', *gamma,'G', *x,'q', lambda));
  double deltaX2(fitRals("Gqr", *gamma,'G', *x,'q', *x,'r', lambda));

  int bcLens[] = { x->lens[0], x->lens[1], x->lens[1] };
  int bcSyms[] = { NS, NS, NS };
  Tensor<complex> bc(3, bcLens, bcSyms, *chi->wrld, "bcRjk", chi->profile);
  bc["Rqr"] = (*x)["Rq"] * (*x)["Rr"];
  (*chi0)["Gqr"] = (*gamma)["RG"] * bc["Rqr"];
  (*chi0)["Gqr"] -= (*chi)["Gqr"];
  double R(frobeniusNorm(*chi0));
  LOG(0) << "R(XXG-chi)=" << R*R << std::endl;
  LOG(3) << "R(X1'-X1)=" << deltaX1 << std::endl;
  LOG(3) << "R(X2'-X2)=" << deltaX2 << std::endl;
  LOG(3) << "R(G'-G)=" << deltaG << std::endl;
  (*chi0)["Gqr"] += (*chi)["Gqr"];
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
  LOG(4) << "inverting Gramian..." << std::endl;
  IterativePseudoInverse<complex> gramianInverse(gramian);
  int bcLens[] = { b.lens[0], b.lens[1], c.lens[1] };
  int bcSyms[] = { NS, NS, NS };
  Tensor<complex> bc(3, bcLens, bcSyms, *chi->wrld, "bcRjk", chi->profile);
  bc["Sjk"] = b["Sj"] * c["Sk"];
  bc["Rjk"] = bc["Sjk"] * gramianInverse.get()["SR"];
  char const indicesA[] = { 'R', idxA, 0 };
  char const indicesBC[] = { 'R', idxB, idxC, 0 };
  LOG(4) << "building new estimate: A[" << indicesA <<
    "] = chi[" << indicesChi << "] * bc[" << indicesBC << "] ..." << std::endl;
  a.contract(1.0, *chi,indicesChi, bc,indicesBC, 0.0, indicesA, fDot);
  LOG(4) << "done" << std::endl;
}

double RalsFtodRankDecomposition::fitRals(
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
  LOG(4) << "inverting Gramian..." << std::endl;
  IterativePseudoInverse<complex> gramianInverse(gramian);
  Tensor<complex> oldA(a);

  int bcLens[] = { b.lens[0], b.lens[1], c.lens[1] };
  int bcSyms[] = { NS, NS, NS };
  LOG(4) << "building outer product..." << std::endl;
  Tensor<complex> bc(3, bcLens, bcSyms, *chi->wrld, "bcRjk", chi->profile);
//  bc["Sjk"] = b["Sj"] * c["Sk"];

  char const indicesA[] = { 'R', idxA, 0 };
  char const indicesBC[] = { 'R', idxB, idxC, 0 };
  char const indicesB[] = { 'R', idxB, 0 };
  char const indicesC[] = { 'R', idxC, 0 };
  LOG(4) << "allocating aux matrices..." << std::endl;
  // a = lambda*a + chi * conj(b*c)
  LOG(4) << "allocating conjB" << std::endl;
  Tensor<complex> conjB(b);
  LOG(4) << "allocating conjC" << std::endl;
  Tensor<complex> conjC(c);
  LOG(4) << "allocating fConj" << std::endl;
  Univar_Function<complex> fConj(&conj<complex>);
// FIXME> cleanup
  // conjB["Rj"] = 0*conjB["Rj"] + conj(b["Rj"])
  LOG(4) << "conjugating B[" << indicesB << "]" << std::endl;
  conjB.sum(1.0, b,"Rj", 0.0,"Rj", fConj); 
  LOG(4) << "conjugating C[" << indicesC << "]" << std::endl;
  conjC.sum(1.0, c,"Rk", 0.0,"Rk", fConj);
//  a.contract(1.0, *chi,indicesChi, bc,indicesBC, lambda, indicesA, fDot);
  LOG(4) << "applying outer product..." << std::endl;
  LOG(4) << "a[" << indicesA << "] = T[" << indicesChi << "] * b * c" << std::endl;
  bc["Sjk"] = conjB["Sj"] * conjC["Sk"];
// FIXME: add to DONT'S
//  a[indicesA] = (*chi)[indicesChi] * (conjB[indicesB] * conjC[indicesC]);
  a[indicesA] *= lambda;
  a[indicesA] += (*chi)[indicesChi] * bc[indicesBC];
  LOG(4) << "applying inverse of Gramian..." << std::endl;
  a.contract(1.0, a,"Si", gramianInverse.get(), "SR", 0.0, "Ri", fDot);
  oldA["Ri"] -= a["Ri"];
  double norm(frobeniusNorm(oldA));
  return norm*norm;
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
  IterativePseudoInverse<complex> gramianInverse(gramian);
  bg["jR"] = b["jS"] * gramianInverse.get()["SR"];
  a.contract(1.0, t,"ij", bg,"jR", 0.0,"iR", fDot);

  Matrix<complex> ab(5,7, NS, *world, "Tij");
  ab["ij"] = a["iR"] * b["jR"];
  ab["ij"] -= t["ij"];
  logMatrix(2, ab);
  double error(frobeniusNorm(ab));
  LOG(2) << "error=" << error << std::endl;
}

