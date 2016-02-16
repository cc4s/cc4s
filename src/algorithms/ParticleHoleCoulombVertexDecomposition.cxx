/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <algorithms/RalsParticleHoleCoulombVertexDecomposition.hpp>
#include <math/CanonicalPolyadicDecomposition.hpp>
// #include <math/ComplexTensor.hpp>
#include <math/RandomTensor.hpp>
#include <math/MathFunctions.hpp>
#include <math/IterativePseudoInverse.hpp>
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
  LOG(0, "RALS") << "Tensor rank decomposition with rank=" << rank << std::endl;
  LOG(1, "RALS") << "decompising Coulomb vertex with NG=" << NG
    << " No=" << No << " Nv=" << Nv << std::endl;

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
  regularizationEstimatorPiiR =
    new AlternatingLeastSquaresRegularizationEstimator(
      swampingThreshold, regularizationFriction, 1
    );
  regularizationEstimatorPiaR =
    new AlternatingLeastSquaresRegularizationEstimator(
      swampingThreshold, regularizationFriction, 1
    );
  regularizationEstimatorLambdaGR =
    new AlternatingLeastSquaresRegularizationEstimator(
      swampingThreshold, regularizationFriction, 1
    );
  int64_t iterationsCount(0);
  int64_t maxIterationsCount(
    getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS)
  );
  double delta(getRealArgument("delta", DEFAULT_DELTA));
  residuum = std::numeric_limits<double>::infinity();
  while (iterationsCount < maxIterationsCount && residuum > delta) {
    fit(iterationsCount);
    ++iterationsCount;
  }
}

void RalsParticleHoleCoulombVertexDecomposition::fit(
  int64_t const iterationsCount
) {
  fitRegularizedAlternatingLeastSquaresFactor(
    *GammaGai,"Gai", *PiaR,'a', *LambdaGR,'G',
    *PiiR,'i', *regularizationEstimatorPiiR
  );
//  realizePi(*PiiR); // normalizePi(*PiiR);
  fitRegularizedAlternatingLeastSquaresFactor(
    *GammaGai,"Gai", *LambdaGR,'G', *PiiR,'i',
    *PiaR,'a', *regularizationEstimatorPiaR
  );
//  realizePi(*PiaR); // normalizePi(*PiaR);
  fitRegularizedAlternatingLeastSquaresFactor(
    *GammaGai,"Gai", *PiiR,'i',*PiaR,'a',
    *LambdaGR,'G', *regularizationEstimatorLambdaGR
  );

  composeCanonicalPolyadicDecompositionTensors(
    *LambdaGR, *PiaR, *PiiR, *Gamma0Gai
  );

  (*Gamma0Gai)["Gai"] -= (*GammaGai)["Gai"];
  residuum = frobeniusNorm(*Gamma0Gai);
  LOG(0, "RALS") << "iteration=" << (iterationsCount+1)
    << " Delta=" << residuum << std::endl;
  (*Gamma0Gai)["Gai"] += (*GammaGai)["Gai"];
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

