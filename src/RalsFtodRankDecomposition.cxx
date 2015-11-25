/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "RalsFtodRankDecomposition.hpp"
#include "Exception.hpp"
#include "util/Log.hpp"
#include "util/ComplexTensor.hpp"
#include "util/RandomTensor.hpp"
#include "util/MathFunctions.hpp"
#include "util/IterativePseudoInverter.hpp"
#include <iostream>
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
  setRandomTensor(*gamma);
  while (true) {
    fit();
  }
}

void RalsFtodRankDecomposition::fit() {
  fitAls("Gqr", *x,'q', *x,'r', *gamma,'G');
  fitAls("Gqr", *x,'r', *gamma,'G', *x,'q');
//  fitAls("Gqr", *gamma,'G', *x,'q', *x,'r');

  int bcLens[] = { x->lens[0], x->lens[1], x->lens[1] };
  int bcSyms[] = { NS, NS, NS };
  Tensor<complex> bc(3, bcLens, bcSyms, *chi->wrld, "bcRjk", chi->profile);
  bc["Rqr"] = (*x)["Rq"] * (*x)["Rr"];
  Tensor<complex> chi0(3, chi->lens, chi->sym, *chi->wrld, "chi0Gqr", chi->profile);
  chi0["Gqr"] = (*gamma)["RG"] * bc["Rqr"];
  chi0["Gqr"] -= (*chi)["Gqr"];
  Scalar<complex> s(*chi->wrld);
  Bivar_Function<complex> fRealDot(&realDot<complex>);
  s.contract(1.0, chi0,"Gqr", chi0,"Gqr", 0.0,"", fRealDot);
  double R(std::real(s.get_val()));
  LOG(4) << R << std::endl;
}

void RalsFtodRankDecomposition::fitAls(
  char const *indicesChi,
  Tensor<complex> &b, char const idxB, Tensor<complex> &c, char const idxC,
  Tensor<complex> &a, char const idxA
) {
  Matrix<complex> gramian(rank, rank, NS, *chi->wrld, "gRS", chi->profile);
  gramian["RS"] =  b["Rj"]*b["Sj"];
  gramian["RS"] *= c["Rk"]*c["Sk"];
  IterativePseudoInverter<complex> gramianInverter(gramian);
  int bcLens[] = { b.lens[0], b.lens[1], c.lens[1] };
  int bcSyms[] = { NS, NS, NS };
  Tensor<complex> bc(3, bcLens, bcSyms, *chi->wrld, "bcRjk", chi->profile);
  bc["Sjk"] = b["Sj"] * c["Sk"];
  bc["Rjk"] = bc["Sjk"] * gramianInverter.invert()["SR"];
  char const indicesA[] = { 'R', idxA, 0 };
  char const indicesBC[] = { 'S', idxB, idxC, 0 };
  a[indicesA] = (*chi)[indicesChi] * bc[indicesBC];
}

