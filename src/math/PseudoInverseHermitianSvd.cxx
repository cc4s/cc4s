#include <math/PseudoInverseHermitianSvd.hpp>

#include <math/MathFunctions.hpp>
#include <util/BlacsWorld.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <util/ScaLapackHermitianEigenSystemDc.hpp>
#include <util/Log.hpp>

using namespace cc4s;
using namespace CTF;


template <typename F>
PseudoInverseHermitianSvd<F>::PseudoInverseHermitianSvd(
  Matrix<F> &A, double epsilon
):
  inverse(A)
{
  // convert CTF matrices into ScaLapack matrices
  BlacsWorld blacsWorld(A.wrld->rank, A.wrld->np);
  ScaLapackMatrix<F> ScaA(A, &blacsWorld);
  ScaLapackMatrix<F> ScaU(ScaA);

  // solve hermitian eigenvectors problem using ScaLapack
  ScaLapackHermitianEigenSystemDc<F> eigenSystem(&ScaA, &ScaU);
  double *lambda(new double[ScaA.lens[0]]);
  eigenSystem.solve(lambda);

  Vector<F> D(A.lens[0], *A.wrld, "Lambda");
  int localLambdaCount(A.wrld->rank == 0 ? A.lens[0] : 0);
  F *lambdaValues(new F[localLambdaCount]);
  int64_t *lambdaIndices(new int64_t[localLambdaCount]);
  bool smallLambda(false);
  for (int64_t i(0); i < localLambdaCount; ++i) {
    lambdaIndices[i] = i;
    smallLambda |= std::abs(lambda[i]) < 1e3*epsilon;
    lambdaValues[i] = (std::abs(lambda[i]) > epsilon) ? 1/lambda[i] : 0;
  }
  if (smallLambda) {
    LOG(1, "PseudoInverseHermitianSvd")
      << "WARNING: Very small eigenvalues occurred" << std::endl;
  }
  D.write(localLambdaCount, lambdaIndices, lambdaValues);
  Matrix<F> U(A);
  ScaU.write(U);
  delete[] lambdaIndices; delete[] lambdaValues;
  delete[] lambda;

  // compose the pseudo inverse from S and the unitary transforms
  inverse["ij"] = U["ik"] * D["k"] * U["kj"];
}

template <typename F>
Matrix<F> &PseudoInverseHermitianSvd<F>::get() {
  return inverse;
}

// instantiate
template
PseudoInverseHermitianSvd<double>::PseudoInverseHermitianSvd(
  Matrix<double> &matrix, double epsilon
);
template
Matrix<double> &PseudoInverseHermitianSvd<double>::get();

template
PseudoInverseHermitianSvd<complex>::PseudoInverseHermitianSvd(
  Matrix<complex> &matrix, double epsilon
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

