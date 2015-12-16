/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <RalsParticleHoleCoulombVertexDecomposition.hpp>
#include <util/ComplexTensor.hpp>
#include <util/RandomTensor.hpp>
#include <util/MathFunctions.hpp>
#include <util/IterativePseudoInverse.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <iostream>
#include <fstream>
#include <random>
#include <limits>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(RalsParticleHoleCoulombVertexDecomposition);

RalsParticleHoleCoulombVertexDecomposition::
  RalsParticleHoleCoulombVertexDecomposition
(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

RalsParticleHoleCoulombVertexDecomposition::
  ~RalsParticleHoleCoulombVertexDecomposition()
{
}

void RalsParticleHoleCoulombVertexDecomposition::run() {
  rank = getIntegerArgument("rank");
  gammaGai = getTensorArgument<complex>("ParticleHoleCoulombVertex");
  double epsilon = getRealArgument("epsilon");
  int nG(gammaGai->lens[0]);
  int nv(gammaGai->lens[1]);
  int no(gammaGai->lens[2]);
  LOG(3) << "rank=" << rank << std::endl;
  LOG(3) << "nG=" << nG << ", no=" << no << ", nv=" << nv << std::endl;
  // allocate factor tensors
  piiR = new Matrix<complex>(
    no, int(rank), NS, *gammaGai->wrld, "PiiR", gammaGai->profile
  );
  piaR = new Matrix<complex>(
    nv, int(rank), NS, *gammaGai->wrld, "PiaR", gammaGai->profile
  );
  lambdaGR = new Matrix<complex>(
    nG, int(rank), NS, *gammaGai->wrld, "LambdaGR", gammaGai->profile
  );
  setRandomTensor(*piiR);
  realizePi(*piiR); normalizePi(*piiR);
  setRandomTensor(*piaR);
  realizePi(*piiR); normalizePi(*piaR);
  setRandomTensor(*lambdaGR);
  allocatedTensorArgument("HoleFactorOrbitals", piiR);
  allocatedTensorArgument("ParticleFactorOrbitals", piaR);
  allocatedTensorArgument("ParticleHoleCoulombFactors", lambdaGR);

  gamma0Gai = new Tensor<complex>(
    3, gammaGai->lens, gammaGai->sym, *gammaGai->wrld, "gamma0Gai",
    gammaGai->profile
  );
  if (isArgumentGiven("ComposedParticleHoleCoulombVertex")) {
    allocatedTensorArgument("ComposedParticleHoleCoulombVertex", gamma0Gai);
  }

  residuum = std::numeric_limits<double>::infinity();
  while (residuum > epsilon) {
    double lambda(getRealArgument("lambda"));
    std::ifstream lambdaFile("lambda");
    if (lambdaFile.is_open()) {
      // override with lambda given from file if existing
      lambdaFile >> lambda;
      lambdaFile.close();
    }
    if (lambda < 0.0) return;
    LOG(2) << "lambda=" << lambda << std::endl;
    fit(lambda);
  }
}

void RalsParticleHoleCoulombVertexDecomposition::fit(double lambda) {
  double deltaPiiR(fitRals("Gai", *piaR,'a', *lambdaGR,'G', *piiR,'i', lambda));
//  realizePi(*piiR); normalizePi(*piiR);
  double deltaPiaR(fitRals("Gai", *lambdaGR,'G', *piiR,'i', *piaR,'a', lambda));
//  realizePi(*piaR); normalizePi(*piaR);
  double deltaLambda(fitRals("Gai", *piiR,'i',*piaR,'a', *lambdaGR,'G',lambda));

  int bcLens[] = { piiR->lens[1], piaR->lens[0], piiR->lens[0] };
  int bcSyms[] = { NS, NS, NS };
  Tensor<complex> bc(
    3, bcLens, bcSyms, *gammaGai->wrld, "bcRjk", gammaGai->profile
  );
  bc["Rai"] = (*piaR)["aR"] * (*piiR)["iR"];
  (*gamma0Gai)["Gai"] = (*lambdaGR)["GR"] * bc["Rai"];
  (*gamma0Gai)["Gai"] -= (*gammaGai)["Gai"];
  residuum = frobeniusNorm(*gamma0Gai);
  residuum *= residuum;
  LOG(0) << "R(Pi_aR Pi_iR Lambda_GR - Gamma^a_iG)=" << residuum << std::endl;
  LOG(3) << "R(Pi_iR'-Pi_iR)=" << deltaPiiR << std::endl;
  LOG(3) << "R(Pi_aR'-Pi_aR)=" << deltaPiaR << std::endl;
  LOG(3) << "R(Lambda_GR'-Lambda_GR)=" << deltaLambda << std::endl;
  (*gamma0Gai)["Gai"] += (*gammaGai)["Gai"];
}

double RalsParticleHoleCoulombVertexDecomposition::fitRals(
  char const *indicesGamma,
  Tensor<complex> &b, char const idxB, Tensor<complex> &c, char const idxC,
  Tensor<complex> &a, char const idxA,
  double lambda
) {
  Matrix<complex> bb(rank, rank, NS, *gammaGai->wrld, "bbRS",gammaGai->profile);
  Matrix<complex> gramian(rank,rank,NS,*gammaGai->wrld,"gRS",gammaGai->profile);
  Bivar_Function<complex> fDot(&dot<complex>);
  LOG(4) << "building Gramian..." << std::endl;
  bb.contract(1.0,b,"jR", b,"jS", 0.0,"SR", fDot);
  gramian.contract(1.0,c,"kR", c,"kS", 0.0,"SR", fDot);
  gramian["SR"] *= bb["SR"];
  gramian["RR"] += lambda;
  LOG(4) << "inverting Gramian..." << std::endl;
  IterativePseudoInverse<complex> gramianInverse(gramian);
  Tensor<complex> oldA(a);

  int bcLens[] = { int(rank), b.lens[0], c.lens[0] };
  int bcSyms[] = { NS, NS, NS };
  LOG(4) << "building outer product..." << std::endl;
  Tensor<complex> bc(3,bcLens,bcSyms,*gammaGai->wrld,"bcRjk",gammaGai->profile);
  char const indicesA[] = { idxA, 'R', 0 };
  char const indicesBC[] = { 'R', idxB, idxC, 0 };
  Tensor<complex> conjB(b);
  Tensor<complex> conjC(c);
  Univar_Function<complex> fConj(&conj<complex>);
  conjB.sum(1.0, b,"jR", 0.0,"jR", fConj); 
  conjC.sum(1.0, c,"kR", 0.0,"kR", fConj);
  bc["Sjk"] = conjB["jS"] * conjC["kS"];
  LOG(4) << "applying outer product..." << std::endl;
  a[indicesA] *= lambda;
  a[indicesA] += (*gammaGai)[indicesGamma] * bc[indicesBC];
  LOG(4) << "applying inverse of Gramian..." << std::endl;
  a.contract(1.0, a,"iS", gramianInverse.get(), "SR", 0.0, "iR", fDot);
  oldA["iR"] -= a["iR"];
  double norm(frobeniusNorm(oldA));
  return norm*norm;
}

/**
 * \brief Normalizes the given factor orbitals.
 */
void RalsParticleHoleCoulombVertexDecomposition::normalizePi(
  Matrix<complex> &pi
) {
  Bivar_Function<complex> fDot(&cc4s::dot<complex>);
  Vector<complex> norm(pi.lens[0], *pi.wrld);
  // norm["q"] = pi["qR"] * conj(pi["qR"])
  norm.contract(1.0, pi,"qR", pi,"qR", 0.0,"q", fDot);
  Matrix<complex> quotient(pi);
  Univar_Function<complex> fSqrt(&cc4s::sqrt<complex>);
  // quotient["qR"] = sqrt(norm["q"])
  quotient.sum(1.0, norm,"q", 0.0,"qR", fSqrt);
  Bivar_Function<complex> fDivide(&cc4s::divide<complex>);
  // pi["qR"] = pi["qR"] / quotient["qR"]
  pi.contract(1.0, pi,"qR", quotient,"qR", 0.0,"qR", fDivide);
}

/**
 * \brief Discards the imaginary part of the given factor orbitals.
 */
void RalsParticleHoleCoulombVertexDecomposition::realizePi(
  Matrix<complex> &pi
) {
  Univar_Function<complex> fConj(&cc4s::conj<complex>);
  Matrix<complex> conjX(pi);
  // conjX["qR"] = pi["qR"]
  conjX.sum(1.0, pi,"qR", 0.0,"qR", fConj);
  pi["qR"] += conjX["qR"];
  pi["qR"] *= 0.5;
}

