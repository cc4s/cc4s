#include <math/PseudoInverseHermitianSvd.hpp>

#include <math/MathFunctions.hpp>
#include <util/BlacsWorld.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <extern/ScaLapack.hpp>
#include <util/Log.hpp>

using namespace cc4s;
using namespace CTF;


template <typename F>
PseudoInverseHermitianSvd<F>::PseudoInverseHermitianSvd(
  Matrix<F> const &matrix_
):
  A(matrix_)
{
}

template <typename F>
Matrix<F> &PseudoInverseHermitianSvd<F>::get() {
  BlacsWorld blacsWorld(A.wrld->rank, A.wrld->np);
  ScaLapackMatrix<F> ScaA(A, &blacsWorld);
  ScaLapackMatrix<F> ScaU(ScaA);
  int offset(1);
  double *lambda(new double[A.lens[0]]);
  int info;
  pheev(
    "EV", "U",
    &ScaA.lens[0],
    ScaA.localValues, &offset, &offset, ScaA.getDescriptor(),
    lambda,
    ScaU.localValues, &offset, &offset, ScaU.getDescriptor(),
    &info
  );
  Vector<F> S(A.lens[0], *A.wrld, "Lambda");
  int localLambdaCount(A.wrld->rank == 0 ? ScaA.lens[0] : 0);
  F *lambdaValues(new F[localLambdaCount]);
  int64_t *lambdaIndices(new int64_t[localLambdaCount]);
  // TODO: invert eigen values
  for (int64_t i(0); i < localLambdaCount; ++i) {
    lambdaIndices[i] = i;
    lambdaValues[i] = lambda[i];
  }
  S.write(localLambdaCount, lambdaIndices, lambdaValues);
  Matrix<F> U(A);
  ScaU.write(U);
  delete[] lambdaIndices; delete[] lambdaValues;
  delete[] lambda;

  Matrix<F> UT(U);
  Univar_Function<F> fConj(conj<F>);
  UT.sum(1.0,U,"ij", 0.0,"ij", fConj);
  U["ij"] = U["ij"] * S["j"];  
  A["ij"] = U["ik"] * UT["jk"];
  return A;
}

// instantiate
template
PseudoInverseHermitianSvd<double>::PseudoInverseHermitianSvd(
  Matrix<double> const &matrix
);
template
Matrix<double> &PseudoInverseHermitianSvd<double>::get();

template
PseudoInverseHermitianSvd<complex>::PseudoInverseHermitianSvd(
  Matrix<complex> const &matrix
);
template
Matrix<complex> &PseudoInverseHermitianSvd<complex>::get();


template <typename F>
DryPseudoInverseHermitianSvd<F>::DryPseudoInverseHermitianSvd(
  DryMatrix<F> const &matrix_
):
  inverse(matrix_)
{
}

template <typename F>
DryMatrix<F> &DryPseudoInverseHermitianSvd<F>::get() {
  return inverse;
}

// instantiate
template
DryPseudoInverseHermitianSvd<double>::DryPseudoInverseHermitianSvd(
  DryMatrix<double> const &matrix
);
template
DryMatrix<double> &DryPseudoInverseHermitianSvd<double>::get();

template
DryPseudoInverseHermitianSvd<complex>::DryPseudoInverseHermitianSvd(
  DryMatrix<complex> const &matrix
);
template
DryMatrix<complex> &DryPseudoInverseHermitianSvd<complex>::get();

