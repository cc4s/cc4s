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
  GammaGai = getTensorArgument<complex>("ParticleHoleCoulombVertex");
  int NG(GammaGai->lens[0]);
  int Nv(GammaGai->lens[1]);
  int No(GammaGai->lens[2]);
  rank = getIntegerArgument("rank", NG);
  LOG(3) << "rank=" << rank << std::endl;
  LOG(3) << "NG=" << NG << ", No=" << No << ", Nv=" << Nv << std::endl;

  // allocate factor tensors
  PiiR = new Matrix<complex>(
    No, int(rank), NS, *GammaGai->wrld, "PiiR", GammaGai->profile
  );
  PiaR = new Matrix<complex>(
    Nv, int(rank), NS, *GammaGai->wrld, "PiaR", GammaGai->profile
  );
  LambdaGR = new Matrix<complex>(
    NG, int(rank), NS, *GammaGai->wrld, "LambdaGR", GammaGai->profile
  );
  setRandomTensor(*PiiR);
  realizePi(*PiiR); normalizePi(*PiiR);
  setRandomTensor(*PiaR);
  realizePi(*PiiR); normalizePi(*PiaR);
  setRandomTensor(*LambdaGR);
  allocatedTensorArgument("HoleFactorOrbitals", PiiR);
  allocatedTensorArgument("ParticleFactorOrbitals", PiaR);
  allocatedTensorArgument("ParticleHoleCoulombFactors", LambdaGR);

  Gamma0Gai = new Tensor<complex>(
    3, GammaGai->lens, GammaGai->sym, *GammaGai->wrld, "Gamma0Gai",
    GammaGai->profile
  );
  if (isArgumentGiven("ComposedParticleHoleCoulombVertex")) {
    allocatedTensorArgument("ComposedParticleHoleCoulombVertex", Gamma0Gai);
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
  LOG(1) << "lambda   s/s_0" << std::endl;
  fitRals(
    "Gai", *PiaR,'a', *LambdaGR,'G', *PiiR,'i', *regularizationEstimatorPiiR
  );
//  realizePi(*PiiR); normalizePi(*PiiR);
  fitRals(
    "Gai", *LambdaGR,'G', *PiiR,'i', *PiaR,'a', *regularizationEstimatorPiaR
  );
//  realizePi(*PiaR); normalizePi(*PiaR);
  fitRals(
    "Gai", *PiiR,'i',*PiaR,'a', *LambdaGR,'G', *regularizationEstimatorLambdaGR
  );

  int bcLens[] = { PiiR->lens[1], PiaR->lens[0], PiiR->lens[0] };
  int bcSyms[] = { NS, NS, NS };
  Tensor<complex> BC(
    3, bcLens, bcSyms, *GammaGai->wrld, "BCRjk", GammaGai->profile
  );
  BC["Rai"] = (*PiaR)["aR"] * (*PiiR)["iR"];
  (*Gamma0Gai)["Gai"] = (*LambdaGR)["GR"] * BC["Rai"];
  (*Gamma0Gai)["Gai"] -= (*GammaGai)["Gai"];
  residuum = frobeniusNorm(*Gamma0Gai);
  LOG(0) << iterationsCount << ":  " << residuum << std::endl;
  (*Gamma0Gai)["Gai"] += (*GammaGai)["Gai"];
}

void RalsParticleHoleCoulombVertexDecomposition::fitRals(
  char const *indicesGamma,
  Tensor<complex> &B, char const idxB, Tensor<complex> &C, char const idxC,
  Tensor<complex> &A, char const idxA,
  RegularizationEstimator &regularizationEstimatorA
) {
  double lambda(regularizationEstimatorA.getLambda());
  Tensor<complex> conjB(B);
  Tensor<complex> conjC(C);
  Univar_Function<complex> fConj(&conj<complex>);
  conjB.sum(1.0, B,"jR", 0.0,"jR", fConj); 
  conjC.sum(1.0, C,"kR", 0.0,"kR", fConj);

  Matrix<complex> BB(rank, rank, NS, *GammaGai->wrld, "BBRS",GammaGai->profile);
  Matrix<complex> gramian(rank,rank,NS,*GammaGai->wrld,"GRS",GammaGai->profile);
  LOG(4) << "building Gramian..." << std::endl;
  BB["SR"] = B["jR"] * conjB["jS"];
  gramian["SR"] = C["kR"] * conjC["kS"];
  gramian["SR"] *= BB["SR"];
  gramian["RR"] += lambda;
  LOG(4) << "inverting Gramian..." << std::endl;
  IterativePseudoInverse<complex> gramianInverse(gramian);

  Tensor<complex> oldA(A);
  char const indicesA[] = { idxA, 'R', 0 };
  char const indicesB[] = { idxB, 'R', 0 };
  char const indicesC[] = { idxC, 'R', 0 };
  LOG(4) << "applying to Gamma..." << std::endl;
  A[indicesA] *= lambda;
  A[indicesA] += (*GammaGai)[indicesGamma] * conjB[indicesB] * conjC[indicesC];
  LOG(4) << "applying inverse of Gramian..." << std::endl;
  Tensor<complex> conjInvGramian(gramianInverse.get());
  conjInvGramian.sum(1.0, conjInvGramian,"SR", 0.0,"SR", fConj);
  A["iR"] = A["iS"] * conjInvGramian["SR"];
  oldA["iR"] -= A["iR"];
  double normDifference(frobeniusNorm(oldA));
  double norm(frobeniusNorm(A));
  double swampingFactor(
    normDifference / norm / regularizationEstimatorA.getSwampingThreshold()
  );
  LOG(1) << lambda << "  " << swampingFactor << std::endl;
  regularizationEstimatorA.update(swampingFactor);
}

/**
 * \brief Normalizes the given factor orbitals.
 */
void RalsParticleHoleCoulombVertexDecomposition::normalizePi(
  Matrix<complex> &Pi
) {
  Bivar_Function<complex> fDot(&cc4s::dot<complex>);
  Vector<complex> norm(Pi.lens[0], *Pi.wrld);
  // norm["q"] = Pi["qR"] * conj(Pi["qR"])
  norm.contract(1.0, Pi,"qR", Pi,"qR", 0.0,"q", fDot);
  Matrix<complex> quotient(Pi);
  Univar_Function<complex> fSqrt(&cc4s::sqrt<complex>);
  // quotient["qR"] = sqrt(norm["q"])
  quotient.sum(1.0, norm,"q", 0.0,"qR", fSqrt);
  Bivar_Function<complex> fDivide(&cc4s::divide<complex>);
  // Pi["qR"] = Pi["qR"] / quotient["qR"]
  Pi.contract(1.0, Pi,"qR", quotient,"qR", 0.0,"qR", fDivide);
}

/**
 * \brief Discards the imaginary part of the given factor orbitals.
 */
void RalsParticleHoleCoulombVertexDecomposition::realizePi(
  Matrix<complex> &Pi
) {
  Univar_Function<complex> fConj(&cc4s::conj<complex>);
  Matrix<complex> conjX(Pi);
  // conjX["qR"] = Pi["qR"]
  conjX.sum(1.0, Pi,"qR", 0.0,"qR", fConj);
  Pi["qR"] += conjX["qR"];
  Pi["qR"] *= 0.5;
}

