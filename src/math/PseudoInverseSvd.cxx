#include <math/PseudoInverseSvd.hpp>

#include <math/MathFunctions.hpp>
#include <util/BlacsWorld.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <util/ScaLapackSingularValueDecomposition.hpp>
#include <util/Log.hpp>

using namespace cc4s;
using namespace CTF;


template <typename F>
PseudoInverseSvd<F>::PseudoInverseSvd(
  Matrix<F> &A, double epsilon
):
  inverse(A)
{
  // convert CTF matrices into ScaLapack matrices
  BlacsWorld blacsWorld(A.wrld->rank, A.wrld->np);
  // TODO: only works for quadratic matrices
  ScaLapackMatrix<F> ScaA(A, &blacsWorld);
  ScaLapackMatrix<F> ScaU(ScaA);
  ScaLapackMatrix<F> ScaVT(ScaA);

  // do SVD using ScaLapack
  ScaLapackSingularValueDecomposition<F> svd(&ScaA, &ScaU, &ScaVT);
  double *sigma(new double[A.lens[0]]);
  svd.decompose(sigma);

  // convert real sigma into complex CTF vector of their pseudo inverses
  // TODO: only works for quadratic matrices
  Vector<F> S(A.lens[0], *A.wrld, "Sigma");
  int localSigmaCount(A.wrld->rank == 0 ? ScaA.lens[0] : 0);
  F *sigmaValues(new F[localSigmaCount]);
  int64_t *sigmaIndices(new int64_t[localSigmaCount]);
  for (int64_t i(0); i < localSigmaCount; ++i) {
    sigmaIndices[i] = i;
    // invert singular values
    sigmaValues[i] = (sigma[i] > epsilon) ? 1 / sigma[i] : 0;
    // TODO: warn about very small singular values above epsilon
  }
  S.write(localSigmaCount, sigmaIndices, sigmaValues);

  // convert ScaLapack result matrices to CTF
  Matrix<F> U(A);
  ScaU.write(U);
  Matrix<F> VT(A);
  ScaVT.write(VT);
  // dump local storage
  delete[] sigmaIndices; delete[] sigmaValues;
  delete[] sigma;

  // recompose in CTF to get pseudo inverse matrix
  inverse["ij"] = VT["ki"] * S["k"] * U["jk"];
  Univar_Function<F> fConj(conj<F>);
  inverse.sum(1.0,inverse,"ij", 0.0,"ij", fConj);
}

template <typename F>
Matrix<F> &PseudoInverseSvd<F>::get() {
  return inverse;
}

// instantiate
template
PseudoInverseSvd<double>::PseudoInverseSvd(
  Matrix<double> &matrix, double epsilon
);
template
Matrix<double> &PseudoInverseSvd<double>::get();

template
PseudoInverseSvd<complex>::PseudoInverseSvd(
  Matrix<complex> &matrix, double epsilon
);
template
Matrix<complex> &PseudoInverseSvd<complex>::get();


template <typename F>
DryPseudoInverseSvd<F>::DryPseudoInverseSvd(
  DryMatrix<F> const &matrix_
):
  inverse(matrix_)
{
}

template <typename F>
DryMatrix<F> &DryPseudoInverseSvd<F>::get() {
  return inverse;
}

// instantiate
template
DryPseudoInverseSvd<double>::DryPseudoInverseSvd(
  DryMatrix<double> const &matrix
);
template
DryMatrix<double> &DryPseudoInverseSvd<double>::get();

template
DryPseudoInverseSvd<complex>::DryPseudoInverseSvd(
  DryMatrix<complex> const &matrix
);
template
DryMatrix<complex> &DryPseudoInverseSvd<complex>::get();

