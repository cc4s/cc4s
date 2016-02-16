/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <algorithms/CoulombVertexDecomposition.hpp>
#include <math/CanonicalPolyadicDecomposition.hpp>
#include <math/RandomTensor.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <limits>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(CoulombVertexDecomposition);

CoulombVertexDecomposition::
  CoulombVertexDecomposition
(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

CoulombVertexDecomposition::
  ~CoulombVertexDecomposition()
{
}

void CoulombVertexDecomposition::run() {
  GammaGqr = getTensorArgument<complex>("CoulombVertex");
  int NG(GammaGqr->lens[0]);
  int Np(GammaGqr->lens[1]);
  rank = getIntegerArgument("rank", NG);
  realFactorOrbitals = getIntegerArgument(
    "realFactorOrbitals", DEFAULT_REAL_FACTOR_ORBITALS
  );
  normalizedFactorOrbitals = getIntegerArgument(
    "normalizedFactorOrbitals", DEFAULT_NORMALIZED_FACTOR_ORBITALS
  );
  LOG(0, "RALS") << "Tensor rank decomposition with rank=" << rank
    << ", realFactorOrbitals=" << realFactorOrbitals
    << ", normalizedFactorOrbitals=" << normalizedFactorOrbitals << std::endl;
  LOG(1, "RALS") << "decompising Coulomb vertex with NG=" << NG
    << " Np=" << Np << std::endl;

  // allocate factor tensors
  PiqR = new Matrix<complex>(
    Np, int(rank), NS, *GammaGqr->wrld, "PiqR", GammaGqr->profile
  );
  LambdaGR = new Matrix<complex>(
    NG, int(rank), NS, *GammaGqr->wrld, "LambdaGR", GammaGqr->profile
  );
  setRandomTensor(*PiqR);
  realizePi(*PiqR); normalizePi(*PiqR);
  setRandomTensor(*LambdaGR);
  allocatedTensorArgument("FactorOrbitals", PiqR);
  allocatedTensorArgument("CoulombFactors", LambdaGR);

  Gamma0Gqr = new Tensor<complex>(
    3, GammaGqr->lens, GammaGqr->sym, *GammaGqr->wrld, "Gamma0Gqr",
    GammaGqr->profile
  );
  if (isArgumentGiven("ComposedCoulombVertex")) {
    allocatedTensorArgument("ComposedCoulombVertex", Gamma0Gqr);
  }

  double swampingThreshold(
    getRealArgument("swampingThreshold", DEFAULT_SWAMPING_THRESHOLD)
  );
  double regularizationFriction(
    getRealArgument("regularizationFriction", DEFAULT_REGULARIZATION_FRICTION)
  );
  regularizationEstimatorPiqR =
    new AlternatingLeastSquaresRegularizationEstimator(
      swampingThreshold, regularizationFriction, 1
    );
  regularizationEstimatorPirR =
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
  Delta = std::numeric_limits<double>::infinity();
  while (iterationsCount < maxIterationsCount && Delta > delta) {
    fit(iterationsCount);
    ++iterationsCount;
  }
}

void CoulombVertexDecomposition::fit(
  int64_t const iterationsCount
) {
  Univar_Function<complex> fConj(&cc4s::conj<complex>);
  fitRegularizedAlternatingLeastSquaresFactor(
    *GammaGqr,"Gqr", *PiqR,'q', *LambdaGR,'G',
    *PirR,'r', *regularizationEstimatorPirR
  );
  if (realFactorOrbitals) realizePi(*PiqR);
  if (normalizedFactorOrbitals) normalizePi(*PiqR);
  // PirR["qR"] = conj(PiqR["qR"])
  PirR->sum(1.0, *PiqR,"qR", 0.0,"qR", fConj);
  fitRegularizedAlternatingLeastSquaresFactor(
    *GammaGqr,"Gqr", *LambdaGR,'G', *PirR,'r',
    *PiqR,'q', *regularizationEstimatorPiqR
  );
  if (realFactorOrbitals) realizePi(*PiqR);
  if (normalizedFactorOrbitals) normalizePi(*PiqR);
  // PiqR["qR"] = conj(PirR["qR"])
  PiqR->sum(1.0, *PirR,"qR", 0.0,"qR", fConj);
  fitRegularizedAlternatingLeastSquaresFactor(
    *GammaGqr,"Gqr", *PirR,'r', *PiqR,'q',
    *LambdaGR,'G', *regularizationEstimatorLambdaGR
  );

  composeCanonicalPolyadicDecompositionTensors(
    *LambdaGR, *PiqR, *PirR, *Gamma0Gqr
  );

  (*Gamma0Gqr)["Gqr"] -= (*GammaGqr)["Gqr"];
  Delta = frobeniusNorm(*Gamma0Gqr);
  LOG(0, "RALS") << "iteration=" << (iterationsCount+1)
    << " Delta=" << Delta << std::endl;
  (*Gamma0Gqr)["Gqr"] += (*GammaGqr)["Gqr"];
}

void CoulombVertexDecomposition::normalizePi(
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

void CoulombVertexDecomposition::realizePi(
  Matrix<complex> &Pi
) {
  Univar_Function<complex> fConj(&cc4s::conj<complex>);
  Matrix<complex> conjX(Pi);
  // conjX["qR"] = conj(Pi["qR"])
  conjX.sum(1.0, Pi,"qR", 0.0,"qR", fConj);
  Pi["qR"] += conjX["qR"];
  Pi["qR"] *= 0.5;
}

