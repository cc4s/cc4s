/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <math/RegularizedAlternatingLeastSquares.hpp>
#include <math/CanonicalPolyadicDecomposition.hpp>
#include <math/MathFunctions.hpp>
#include <math/IterativePseudoInverse.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>

using namespace CTF;
using namespace cc4s;

template <typename F=double>
void cc4s::fitRegularizedAlternatingLeastSquaresFactor(
  Tensor<F> &T, char const *indicesT,
  Tensor<F> &B, char const idxB, Tensor<F> &C, char const idxC,
  Tensor<F> &A, char const idxA,
  AlternatingLeastSquaresRegularizationEstimator &regularizationEstimatorA
) {
  double lambda(regularizationEstimatorA.getLambda());
  Tensor<F> conjB(B);
  Tensor<F> conjC(C);
  Univar_Function<F> fConj(&conj<F>);
  conjB.sum(1.0, B,"jR", 0.0,"jR", fConj); 
  conjC.sum(1.0, C,"kR", 0.0,"kR", fConj);

  Matrix<F> BB(B.lens[1], B.lens[1], NS, *T.wrld, "BBRS", T.profile);
  Matrix<F> gramian(B.lens[1], B.lens[1], NS, *T.wrld,"GRS", T.profile);
  LOG(4, "RALS") << "building Gramian..." << std::endl;
  BB["SR"] = B["jR"] * conjB["jS"];
  gramian["SR"] = C["kR"] * conjC["kS"];
  gramian["SR"] *= BB["SR"];
  gramian["RR"] += lambda;
  LOG(4, "RALS") << "inverting Gramian..." << std::endl;
  IterativePseudoInverse<F> gramianInverse(gramian);

  Tensor<F> previousA(A);
  contractWithCanonicalPolyadicDecompositionTensors(
    T, indicesT, conjB, idxB, conjC, idxC, A, idxA
  );
  A["iR"] += lambda * previousA["iR"];
  LOG(4, "RALS") << "applying inverse of Gramian..." << std::endl;
  Tensor<F> conjInvGramian(gramianInverse.get());
  conjInvGramian.sum(1.0, conjInvGramian,"SR", 0.0,"SR", fConj);
  A["iR"] = A["iS"] * conjInvGramian["SR"];
  previousA["iR"] -= A["iR"];
  double normDifference(frobeniusNorm(previousA));
  double norm(frobeniusNorm(A));
  double swampingFactor(
    normDifference / norm / regularizationEstimatorA.getSwampingThreshold()
  );
  LOG(1, "RALS") << "lambda=" << lambda << " s/s_0=" << swampingFactor
    << std::endl;
  regularizationEstimatorA.update(swampingFactor);
}

// instantiate
template
void cc4s::fitRegularizedAlternatingLeastSquaresFactor(
  Tensor<double> &T, char const *indicesT,
  Tensor<double> &B, char const idxB, Tensor<double> &C, char const idxC,
  Tensor<double> &A, char const idxA,
  AlternatingLeastSquaresRegularizationEstimator &regularizationEstimatorA
);
template
void cc4s::fitRegularizedAlternatingLeastSquaresFactor(
  Tensor<complex> &T, char const *indicesT,
  Tensor<complex> &B, char const idxB, Tensor<complex> &C, char const idxC,
  Tensor<complex> &A, char const idxA,
  AlternatingLeastSquaresRegularizationEstimator &regularizationEstimatorA
);

