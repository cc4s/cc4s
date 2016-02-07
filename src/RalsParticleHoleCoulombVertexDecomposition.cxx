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
  gammaGai = getTensorArgument<complex>("ParticleHoleCoulombVertex");
  int nG(gammaGai->lens[0]);
  int nv(gammaGai->lens[1]);
  int no(gammaGai->lens[2]);
  rank = getIntegerArgument("rank", nG);
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

  double swampingThreshold(
    getRealArgument("swampingThreshold", DEFAULT_SWAMPING_THRESHOLD)
  );
  double regularizationFriction(
    getRealArgument("regularizationFriction", DEFAULT_REGULARIZATION_FRICTION)
  );
  regularizationEstimatorPiiR = new RegularizationEstimator(
    swampingThreshold, regularizationFriction, 1
  );
  regularizationEstimatorPiaR = new RegularizationEstimator(
    swampingThreshold, regularizationFriction, 1
  );
  regularizationEstimatorLambdaGR = new RegularizationEstimator(
    swampingThreshold, regularizationFriction, 1
  );
  int64_t iterationsCount(0);
  int64_t maxIterationsCount(
    getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS)
  );
  double delta(getRealArgument("delta", DEFAULT_DELTA));
  residuum = std::numeric_limits<double>::infinity();
  LOG(0) << "n :  Delta" << std::endl;
  while (iterationsCount < maxIterationsCount && residuum > delta) {
    fit(iterationsCount);
    ++iterationsCount;
  }
}

void RalsParticleHoleCoulombVertexDecomposition::fit(
  int64_t const iterationsCount
) {
  LOG(1) << "lambda   s" << std::endl;
  fitRals(
    "Gai", *piaR,'a', *lambdaGR,'G', *piiR,'i', *regularizationEstimatorPiiR
  );
//  realizePi(*piiR); normalizePi(*piiR);
  fitRals(
    "Gai", *lambdaGR,'G', *piiR,'i', *piaR,'a', *regularizationEstimatorPiaR
  );
//  realizePi(*piaR); normalizePi(*piaR);
  fitRals(
    "Gai", *piiR,'i',*piaR,'a', *lambdaGR,'G', *regularizationEstimatorLambdaGR
  );

  int bcLens[] = { piiR->lens[1], piaR->lens[0], piiR->lens[0] };
  int bcSyms[] = { NS, NS, NS };
  Tensor<complex> bc(
    3, bcLens, bcSyms, *gammaGai->wrld, "bcRjk", gammaGai->profile
  );
  bc["Rai"] = (*piaR)["aR"] * (*piiR)["iR"];
  (*gamma0Gai)["Gai"] = (*lambdaGR)["GR"] * bc["Rai"];
  (*gamma0Gai)["Gai"] -= (*gammaGai)["Gai"];
  residuum = frobeniusNorm(*gamma0Gai);
  LOG(0) << iterationsCount << ":  " << residuum << std::endl;
  (*gamma0Gai)["Gai"] += (*gammaGai)["Gai"];
}

void RalsParticleHoleCoulombVertexDecomposition::fitRals(
  char const *indicesGamma,
  Tensor<complex> &b, char const idxB, Tensor<complex> &c, char const idxC,
  Tensor<complex> &a, char const idxA,
  RegularizationEstimator &regularizationEstimatorA
) {
  double lambda(regularizationEstimatorA.getLambda());
  Tensor<complex> conjB(b);
  Tensor<complex> conjC(c);
  Univar_Function<complex> fConj(&conj<complex>);
  conjB.sum(1.0, b,"jR", 0.0,"jR", fConj); 
  conjC.sum(1.0, c,"kR", 0.0,"kR", fConj);

  Matrix<complex> bb(rank, rank, NS, *gammaGai->wrld, "bbRS",gammaGai->profile);
  Matrix<complex> gramian(rank,rank,NS,*gammaGai->wrld,"gRS",gammaGai->profile);
  LOG(4) << "building Gramian..." << std::endl;
  bb["SR"] = b["jR"] * conjB["jS"];
  gramian["SR"] = c["kR"] * conjC["kS"];
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
  bc["Sjk"] = conjB["jS"] * conjC["kS"];
  LOG(4) << "applying outer product..." << std::endl;
  a[indicesA] *= lambda;
  a[indicesA] += (*gammaGai)[indicesGamma] * bc[indicesBC];
  LOG(4) << "applying inverse of Gramian..." << std::endl;
  Tensor<complex> conjInvGramian(gramianInverse.get());
  conjInvGramian.sum(1.0, conjInvGramian,"SR", 0.0,"SR", fConj);
  a["iR"] = a["iS"] * conjInvGramian["SR"];
  oldA["iR"] -= a["iR"];
  double normDifference(frobeniusNorm(oldA));
  double norm(frobeniusNorm(a));
  double swampingFactor(normDifference / norm);
  LOG(1) << lambda << "  " << swampingFactor << std::endl;
  regularizationEstimatorA.update(swampingFactor);
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

