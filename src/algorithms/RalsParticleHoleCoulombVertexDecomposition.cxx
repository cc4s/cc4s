/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <algorithms/RalsParticleHoleCoulombVertexDecomposition.hpp>
#include <math/ComplexTensor.hpp>
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
  while (iterationsCount < maxIterationsCount && residuum > delta) {
    fit(iterationsCount);
    ++iterationsCount;
  }
}

void RalsParticleHoleCoulombVertexDecomposition::fit(
  int64_t const iterationsCount
) {
  fitRals(
    "Gai", *PiaR,'a', *LambdaGR,'G', *PiiR,'i', *regularizationEstimatorPiiR
  );
//  realizePi(*PiiR); // normalizePi(*PiiR);
  fitRals(
    "Gai", *LambdaGR,'G', *PiiR,'i', *PiaR,'a', *regularizationEstimatorPiaR
  );
//  realizePi(*PiaR); // normalizePi(*PiaR);
  fitRals(
    "Gai", *PiiR,'i',*PiaR,'a', *LambdaGR,'G', *regularizationEstimatorLambdaGR
  );

  calculateGamma0();

  (*Gamma0Gai)["Gai"] -= (*GammaGai)["Gai"];
  residuum = frobeniusNorm(*Gamma0Gai);
  LOG(0, "RALS") << "iteration=" << iterationsCount
    << " Delta=" << residuum << std::endl;
  (*Gamma0Gai)["Gai"] += (*GammaGai)["Gai"];
}


/**
 * \brief Calculates \f$\Lambda^a_{iG}\f$.
 */
void RalsParticleHoleCoulombVertexDecomposition::calculateGamma0() {
  if (PiaR->lens[0] < LambdaGR->lens[0]) {
    int lens[] = { PiiR->lens[1], PiaR->lens[0], PiiR->lens[0] };
    int syms[] = { NS, NS, NS };
    Tensor<complex> PiPi(
      3, lens, syms, *GammaGai->wrld, "PiPiRjk", GammaGai->profile
    );
    PiPi["Rai"] = (*PiaR)["aR"] * (*PiiR)["iR"];
    (*Gamma0Gai)["Gai"] = (*LambdaGR)["GR"] * PiPi["Rai"];
  } else {
    int lens[] = { PiiR->lens[1], LambdaGR->lens[0], PiiR->lens[0] };
    int syms[] = { NS, NS, NS };
    Tensor<complex> LambdaPi(
      3, lens, syms, *GammaGai->wrld, "BCRjk", GammaGai->profile
    );
    LambdaPi["RGi"] = (*LambdaGR)["GR"] * (*PiiR)["iR"];
    (*Gamma0Gai)["Gai"] = (*PiaR)["aR"] * LambdaPi["RGi"];
  }
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
  LOG(4, "RALS") << "building Gramian..." << std::endl;
  BB["SR"] = B["jR"] * conjB["jS"];
  gramian["SR"] = C["kR"] * conjC["kS"];
  gramian["SR"] *= BB["SR"];
  gramian["RR"] += lambda;
  LOG(4, "RALS") << "inverting Gramian..." << std::endl;
  IterativePseudoInverse<complex> gramianInverse(gramian);

  Tensor<complex> oldA(A);
  A["iR"] *= lambda;
  applyToGamma(indicesGamma, conjB, idxB, conjC, idxC, A, idxA);
//  applyToGammaSliced(indicesGamma, conjB, idxB, conjC, idxC, A, idxA);
  LOG(4, "RALS") << "applying inverse of Gramian..." << std::endl;
  Tensor<complex> conjInvGramian(gramianInverse.get());
  conjInvGramian.sum(1.0, conjInvGramian,"SR", 0.0,"SR", fConj);
  A["iR"] = A["iS"] * conjInvGramian["SR"];
  oldA["iR"] -= A["iR"];
  double normDifference(frobeniusNorm(oldA));
  double norm(frobeniusNorm(A));
  double swampingFactor(
    normDifference / norm / regularizationEstimatorA.getSwampingThreshold()
  );
  LOG(1, "RALS") << "lambda=" << lambda << " s/s_0=" << swampingFactor
    << std::endl;
  regularizationEstimatorA.update(swampingFactor);
}


/**
 * \brief Calculates \f$A_{iR} = A_{iR} + T_{ijk}{B^\ast}^{jR}{C^\ast}^{kR}\f$
 * using a contraction order with minimal memory footprint.
 */
void RalsParticleHoleCoulombVertexDecomposition::applyToGamma(
  char const *indicesGamma,
  Tensor<complex> &conjB, char const idxB,
  Tensor<complex> &conjC, char const idxC,
  Tensor<complex> &A, char const idxA
) {
  char const indicesA[] = { idxA, 'R', 0 };
  char const indicesB[] = { idxB, 'R', 0 };
  char const indicesC[] = { idxC, 'R', 0 };
  // choose contraction order with minimal memory footprint
  int largestIndex(
    std::max(
      std::max(GammaGai->lens[0], GammaGai->lens[1]), GammaGai->lens[2]
    )
  );
  if (A.lens[0] == largestIndex) {
    // A has largest index: contract conjB and conjC first
    LOG(4, "RALS") << "applying to Gamma with largest A..." << std::endl;
    const char indicesBC[] = { idxB, idxC, 'R' , 0};
    int lens[] = { conjB.lens[0], conjC.lens[0], A.lens[1] };
    int syms[] = { NS, NS, NS };
    Tensor<complex> BC(3, lens, syms, *GammaGai->wrld, "BC");
    BC[indicesBC] = conjB[indicesB]*conjC[indicesC];
    A[indicesA] += (*GammaGai)[indicesGamma] * BC[indicesBC];
  } else if (conjB.lens[0] == largestIndex) {
    // B has largest index: contract Gamma and conjB first
    LOG(4, "RALS") << "applying to Gamma with largest B..." << std::endl;
    const char indicesGammaB[] = { idxA, idxC, 'R' , 0};
    int lens[] = { A.lens[0], conjC.lens[0], A.lens[1] };
    int syms[] = { NS, NS, NS };
    Tensor<complex> GammaB(3, lens, syms, *GammaGai->wrld, "GammaB");
    GammaB[indicesGammaB] = (*GammaGai)[indicesGamma]*conjB[indicesB];
    A[indicesA] += GammaB[indicesGammaB] * conjC[indicesC];
  } else {
    // C has largest index: contract Gamma and conjC first
    LOG(4, "RALS") << "applying to Gamma with largest C..." << std::endl;
    const char indicesGammaC[] = { idxA, idxB, 'R' , 0};
    int lens[] = { A.lens[0], conjB.lens[0], A.lens[1] };
    int syms[] = { NS, NS, NS };
    Tensor<complex> GammaC(3, lens, syms, *GammaGai->wrld, "GammaC");
    GammaC[indicesGammaC] = (*GammaGai)[indicesGamma]*conjC[indicesC];
    A[indicesA] += GammaC[indicesGammaC] * conjB[indicesB];
  }
}

void RalsParticleHoleCoulombVertexDecomposition::applyToGammaSliced(
  char const *indicesGamma,
  Tensor<complex> &conjB, char const idxB,
  Tensor<complex> &conjC, char const idxC,
  Tensor<complex> &A, char const idxA
) {
  char const indicesA[] = { idxA, 'S', 0 };
  char const indicesB[] = { idxB, 'S', 0 };
  char const indicesC[] = { idxC, 'S', 0 };
  int dimA(std::string(indicesGamma).find(idxA));
  int dimB(std::string(indicesGamma).find(idxB));
  int dimC(std::string(indicesGamma).find(idxC));
  int nA(GammaGai->lens[dimA]);
  int nB(GammaGai->lens[dimB]), nC(GammaGai->lens[dimC]);
  int nR(conjB.lens[1]);
  int contractionWindow(getIntegerArgument("contractionWindow", 32));
  for (int i(0); i < nA; i += contractionWindow) {
//    for (int j(0); j < nB; j += contractionWindow) {
      int GammaGaiStart[3], GammaGaiEnd[3];
      GammaGaiStart[dimA] = i;
      GammaGaiEnd[dimA] = std::min(i+contractionWindow, nA);
//      GammaGaiStart[dimB] = j;
//      GammaGaiEnd[dimB] = std::min(j+contractionWindow, nB);
      GammaGaiStart[dimB] = 0;
      GammaGaiEnd[dimB] = nB;
      GammaGaiStart[dimC] = 0;
      GammaGaiEnd[dimC] = nC;
      int TCLens[] = {
        GammaGaiEnd[dimA] - GammaGaiStart[dimA],
        GammaGaiEnd[dimB] - GammaGaiStart[dimB],
        nR
      };
      int TCSyms[] = { NS, NS, NS };
      Tensor<complex> TC(3, TCLens, TCSyms, *conjB.wrld, "TCijS");
      char const indicesTC[] = { idxA, idxB, 'S', 0 };
      LOG(4, "RALS") << "slicing Gamma..." << std::endl;
      TC[indicesTC] =
        GammaGai->slice(GammaGaiStart,GammaGaiEnd)[indicesGamma] *
        conjC[indicesC];
//      int BStart[] = { GammaGaiStart[dimB], 0 };
//      int BEnd[] = { GammaGaiEnd[dimB], nR };
      int TBCLens[] = { GammaGaiEnd[dimA]-GammaGaiStart[dimA], nR };
      int TBCSyms[] = { NS, NS };
      Tensor<complex> TBC(2, TBCLens, TBCSyms, *conjB.wrld, "TBCijS");
//      LOG(4, "RALS") << "slicing B..." << std::endl;
//      TBC[indicesA] = TC[indicesTC] * conjB.slice(BStart,BEnd)[indicesB];
      TBC[indicesA] = TC[indicesTC] * conjB[indicesB];
      int AStart[] = { GammaGaiStart[dimA], 0 };
      int AEnd[] = { GammaGaiEnd[dimA], nR };
      int TBCOffsets[] = { 0, 0, 0 };
      LOG(4, "RALS") << "slicing into A..." << std::endl;
      A.slice(AStart,AEnd,1.0, TBC,TBCOffsets,TBCLens,1.0);
//    }
  }
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

