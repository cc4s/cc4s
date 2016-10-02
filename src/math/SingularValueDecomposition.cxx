#include <math/SingularValueDecomposition.hpp>

#include <math/MathFunctions.hpp>
#include <util/BlacsWorld.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <util/Log.hpp>

using namespace cc4s;
using namespace CTF;


template <typename F>
SingularValueDecomposition<F>::SingularValueDecomposition(
  Matrix<F> const &matrix_
):
  inverse(matrix_)
{
}

template <typename F>
Matrix<F> &SingularValueDecomposition<F>::get(){
  BlacsWorld blacsWorld(inverse.wrld->rank, inverse.wrld->np);
  ScaLapackMatrix<F> scaLapackA(inverse, &blacsWorld);
/*
  int n(inverse.lens[0]);
  int iA, jA;
  int descA[9];
  complex *a;
  double *sigma;
  complex *U, *VT;
  int iU, jU, iVT, jVT;
  int descU[9], descVT[9];
  complex *work;
  int workCount;
  double *realWork;
  int info;
  pzgesvd_(
    "V", "V", &n, &n, a, &iA, &jA, descA,
    sigma,
    U, &iU, &jU, descU,
    VT, &iVT, &jVT, descVT,
    work, &workCount, realWork,
    &info
  );
*/
  return inverse;
}

// instantiate
template
SingularValueDecomposition<double>::SingularValueDecomposition(
  Matrix<double> const &matrix
);
template
Matrix<double> &SingularValueDecomposition<double>::get();

template
SingularValueDecomposition<complex>::SingularValueDecomposition(
  Matrix<complex> const &matrix
);
template
Matrix<complex> &SingularValueDecomposition<complex>::get();


template <typename F>
DrySingularValueDecomposition<F>::DrySingularValueDecomposition(
  DryMatrix<F> const &matrix_
):
  inverse(matrix_)
{
}

template <typename F>
DryMatrix<F> &DrySingularValueDecomposition<F>::get() {
  return inverse;
}

// instantiate
template
DrySingularValueDecomposition<double>::DrySingularValueDecomposition(
  DryMatrix<double> const &matrix
);
template
DryMatrix<double> &DrySingularValueDecomposition<double>::get();

template
DrySingularValueDecomposition<complex>::DrySingularValueDecomposition(
  DryMatrix<complex> const &matrix
);
template
DryMatrix<complex> &DrySingularValueDecomposition<complex>::get();

