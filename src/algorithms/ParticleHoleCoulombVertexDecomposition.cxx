/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <algorithms/ParticleHoleCoulombVertexDecomposition.hpp>
#include <math/CanonicalPolyadicDecomposition.hpp>
#include <math/RandomTensor.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
// FIXME: only used by MP2 functions which need to be factored.
#include <math/ComplexTensor.hpp>
#include <limits>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(ParticleHoleCoulombVertexDecomposition);

ParticleHoleCoulombVertexDecomposition::
  ParticleHoleCoulombVertexDecomposition
(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ParticleHoleCoulombVertexDecomposition::
  ~ParticleHoleCoulombVertexDecomposition()
{
}

void ParticleHoleCoulombVertexDecomposition::run() {
  GammaGai = getTensorArgument<complex>("ParticleHoleCoulombVertex");
  int NG(GammaGai->lens[0]);
  int Nv(GammaGai->lens[1]);
  int No(GammaGai->lens[2]);
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
  LOG(1, "RALS") << "decomposing Coulomb vertex with NG=" << NG
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
  allocatedTensorArgument<complex>("HoleFactorOrbitals", PiiR);
  allocatedTensorArgument<complex>("ParticleFactorOrbitals", PiaR);
  allocatedTensorArgument<complex>("ParticleHoleCoulombFactors", LambdaGR);

  Gamma0Gai = new Tensor<complex>(
    3, GammaGai->lens, GammaGai->sym, *GammaGai->wrld, "Gamma0Gai",
    GammaGai->profile
  );
  if (isArgumentGiven("ComposedParticleHoleCoulombVertex")) {
    allocatedTensorArgument<complex>(
      "ComposedParticleHoleCoulombVertex", Gamma0Gai
    );
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
  Delta = std::numeric_limits<double>::infinity();
  epsilonStep = getIntegerArgument("epsilonStep", DEFAULT_EPSILON_STEP);
  if (epsilonStep > 0) mp2Energy = evaluateMp2(*GammaGai);
  while (iterationsCount < maxIterationsCount && Delta > delta) {
    fit(iterationsCount);
    ++iterationsCount;
    if (epsilonStep > 0 && iterationsCount%epsilonStep == 0) evaluateMp2Error();
  }
}

void ParticleHoleCoulombVertexDecomposition::dryRun() {
  // NOTE that in the dry run GammaGai,... are local variables
  DryTensor<complex> *GammaGai(
    getTensorArgument<complex, DryTensor<complex>>("ParticleHoleCoulombVertex")
  );
  int NG(GammaGai->lens[0]);
  int Nv(GammaGai->lens[1]);
  int No(GammaGai->lens[2]);
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
    << " No=" << No << " Nv=" << Nv << std::endl;

  // allocate factor tensors
  DryTensor<complex> *PiiR(new DryMatrix<complex>(No, int(rank), NS));
  DryTensor<complex> *PiaR(new DryMatrix<complex>(Nv, int(rank), NS));
  DryTensor<complex> *LambdaGR(new DryMatrix<complex>(NG, int(rank), NS));
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "HoleFactorOrbitals", PiiR
  );
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "ParticleFactorOrbitals", PiaR
  );
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "ParticleHoleCoulombFactors", LambdaGR
  );

  DryTensor<complex> *Gamma0Gai(GammaGai);
  if (isArgumentGiven("ComposedParticleHoleCoulombVertex")) {
    allocatedTensorArgument<complex, DryTensor<complex>>(
      "ComposedParticleHoleCoulombVertex", Gamma0Gai
    );
  }

  dryFit(GammaGai, PiaR, PiiR, LambdaGR, Gamma0Gai);

  epsilonStep = getIntegerArgument("epsilonStep", DEFAULT_EPSILON_STEP);
  if (epsilonStep > 0) dryEvaluateMp2(*GammaGai);
}


void ParticleHoleCoulombVertexDecomposition::fit(
  int64_t const iterationsCount
) {
  fitRegularizedAlternatingLeastSquaresFactor(
    *GammaGai,"Gai", *PiaR,'a', *LambdaGR,'G',
    *PiiR,'i', regularizationEstimatorPiiR
  );
  if (realFactorOrbitals) realizePi(*PiiR);
  if (normalizedFactorOrbitals) normalizePi(*PiiR);
  fitRegularizedAlternatingLeastSquaresFactor(
    *GammaGai,"Gai", *LambdaGR,'G', *PiiR,'i',
    *PiaR,'a', regularizationEstimatorPiaR
  );
  if (realFactorOrbitals) realizePi(*PiaR);
  if (normalizedFactorOrbitals) normalizePi(*PiaR);
  fitRegularizedAlternatingLeastSquaresFactor(
    *GammaGai,"Gai", *PiiR,'i',*PiaR,'a',
    *LambdaGR,'G', regularizationEstimatorLambdaGR
  );

  composeCanonicalPolyadicDecompositionTensors(
    *LambdaGR, *PiaR, *PiiR, *Gamma0Gai
  );

  (*Gamma0Gai)["Gai"] -= (*GammaGai)["Gai"];
  Delta = frobeniusNorm(*Gamma0Gai);
  LOG(0, "RALS") << "iteration=" << (iterationsCount+1)
    << " Delta=" << Delta << std::endl;
  (*Gamma0Gai)["Gai"] += (*GammaGai)["Gai"];
}

void ParticleHoleCoulombVertexDecomposition::dryFit(
  DryTensor<complex> *GammaGai,
  DryTensor<complex> *PiaR,
  DryTensor<complex> *PiiR,
  DryTensor<complex> *LambdaGR,
  DryTensor<complex> *Gamma0Gai
) {
  dryFitRegularizedAlternatingLeastSquaresFactor(
    *GammaGai,"Gai", *PiaR,'a', *LambdaGR,'G',
    *PiiR,'i'
  );
  dryFitRegularizedAlternatingLeastSquaresFactor(
    *GammaGai,"Gai", *LambdaGR,'G', *PiiR,'i',
    *PiaR,'a'
  );
  dryFitRegularizedAlternatingLeastSquaresFactor(
    *GammaGai,"Gai", *PiiR,'i',*PiaR,'a',
    *LambdaGR,'G'
  );

  dryComposeCanonicalPolyadicDecompositionTensors(
    *LambdaGR, *PiaR, *PiiR, *Gamma0Gai
  );
}


void ParticleHoleCoulombVertexDecomposition::normalizePi(
  Tensor<complex> &Pi
) {
  Bivar_Function<complex> fDot(&cc4s::dot<complex>);
  Vector<complex> norm(Pi.lens[0], *Pi.wrld);
  // norm["q"] = Pi["qR"] * conj(Pi["qR"])
  norm.contract(1.0, Pi,"qR", Pi,"qR", 0.0,"q", fDot);
  Tensor<complex> quotient(Pi);
  Univar_Function<complex> fSqrt(&cc4s::sqrt<complex>);
  // quotient["qR"] = sqrt(norm["q"])
  quotient.sum(1.0, norm,"q", 0.0,"qR", fSqrt);
  Bivar_Function<complex> fDivide(&cc4s::divide<complex>);
  // Pi["qR"] = Pi["qR"] / quotient["qR"]
  Pi.contract(1.0, Pi,"qR", quotient,"qR", 0.0,"qR", fDivide);
  double normalization(frobeniusNorm(quotient));
  LOG(2, "RALS") << "|normalization quotient|=" << normalization << std::endl;
}

void ParticleHoleCoulombVertexDecomposition::realizePi(
  Tensor<complex> &Pi
) {
  Univar_Function<complex> fConj(&cc4s::conj<complex>);
  Tensor<complex> conjX(Pi);
  // conjX["qR"] = Pi["qR"]
  conjX.sum(1.0, Pi,"qR", 0.0,"qR", fConj);
  Pi["qR"] += conjX["qR"];
  Pi["qR"] *= 0.5;
}

void ParticleHoleCoulombVertexDecomposition::evaluateMp2Error() {
  double approximateMp2(evaluateMp2(*Gamma0Gai));
  LOG(2, "RALS") << "approxiomate mp2=" << approximateMp2 << std::endl;
  LOG(2, "RALS") << "mp2=" << mp2Energy << std::endl;
  LOG(1, "RALS") << "epsilon=" << ((approximateMp2-mp2Energy)/mp2Energy) << std::endl;
}

double ParticleHoleCoulombVertexDecomposition::evaluateMp2(
  Tensor<complex> &Gamma
) {
  // FIXME: use sequence type for output to factorize this code
  Tensor<> realGammaGai(
    3, Gamma.lens, Gamma.sym, *Gamma.wrld, "RealGammaGai"
  );
  Tensor<> imagGammaGai(
    3, Gamma.lens, Gamma.sym, *Gamma.wrld, "ImagGammaGai"
  );
  // split into real and imaginary parts
  fromComplexTensor(Gamma, realGammaGai, imagGammaGai);

  // allocate coulomb integrals
  int Nv(Gamma.lens[1]);
  int No(Gamma.lens[2]);
  int lens[] = { Nv, Nv, No, No };
  int syms[] = { NS, NS, NS, NS };
  Tensor<> Vabij(4, lens, syms, *Gamma.wrld, "Vabij");
  Vabij["abij"] =  realGammaGai["gai"] * realGammaGai["gbj"];
  Vabij["abij"] += imagGammaGai["gai"] * imagGammaGai["gbj"];

  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
 
  Tensor<> Tabij(false, Vabij);
  Tabij["abij"] =  (*epsi)["i"];
  Tabij["abij"] += (*epsi)["j"];
  Tabij["abij"] -= (*epsa)["a"];
  Tabij["abij"] -= (*epsa)["b"];

  Bivar_Function<> fDivide(&divide<double>);
  Tabij.contract(1.0, Vabij,"abij", Tabij,"abij", 0.0,"abij", fDivide);

  Scalar<> energy(*Vabij.wrld);
  double dire, exce;

  energy[""] = 2.0 * Tabij["abij"] * Vabij["abij"];
  dire = energy.get_val();
  energy[""] = Tabij["abji"] * Vabij["abij"];
  exce = -1.0 * energy.get_val();
  return dire + exce;
}


void ParticleHoleCoulombVertexDecomposition::dryEvaluateMp2(
  DryTensor<complex> &Gamma
) {
  // FIXME: use sequence type for output to factorize this code
  DryTensor<> realGammaGai(3, Gamma.lens.data(), Gamma.syms.data());
  DryTensor<> imagGammaGai(3, Gamma.lens.data(), Gamma.syms.data());

  // allocate coulomb integrals
  int Nv(Gamma.lens[1]);
  int No(Gamma.lens[2]);
  int lens[] = { Nv, Nv, No, No };
  int syms[] = { NS, NS, NS, NS };
  DryTensor<> Vabij(4, lens, syms);

  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies")
  );
  epsi->use();
  epsa->use();
  DryTensor<> Tabij(Vabij);
  DryScalar<> energy();
}

