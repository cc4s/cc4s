#include <math/SingularValueDecomposition.hpp>

#include <math/MathFunctions.hpp>
#include <util/BlacsWorld.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <extern/ScaLapack.hpp>
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
Matrix<F> &SingularValueDecomposition<F>::get() {
  BlacsWorld blacsWorld(inverse.wrld->rank, inverse.wrld->np);
  Matrix<F> I(inverse);
  ScaLapackMatrix<F> A(inverse, &blacsWorld);
  ScaLapackMatrix<F> U(A);
  ScaLapackMatrix<F> VT(A);
  int offset(1);
/*
  pgemm(
    "None", "None",
    &A.lens[0], &A.lens[1], &A.lens[1],
    &alpha,
    A.localValues, &offset, &offset, A.getDescriptor(),
    A.localValues, &offset, &offset, A.getDescriptor(),
    &beta,
    C.localValues, &offset, &offset, C.getDescriptor()
  );
*/
  double *sigma(new double[A.lens[0]]);
  F optimalWork;
  double optimalRealWork;
  int info;
  int workCount(-1);
  pgesvd(
    "V", "V",
    &A.lens[0], &A.lens[1],
    A.localValues, &offset, &offset, A.getDescriptor(),
    sigma,
    U.localValues, &offset, &offset, U.getDescriptor(),
    VT.localValues, &offset, &offset, VT.getDescriptor(),
    &optimalWork, &workCount, &optimalRealWork,
    &info
  );
  LOG(4, "SVD") << "work size=" << optimalWork << "," << optimalRealWork << std::endl;
  workCount = static_cast<int>(std::real(optimalWork)+0.5);
  F *work(new F[workCount]);
  double *realWork(new double[static_cast<int64_t>(optimalRealWork+0.5)]);
  pgesvd(
    "V", "V",
    &A.lens[0], &A.lens[1],
    A.localValues, &offset, &offset, A.getDescriptor(),
    sigma,
    U.localValues, &offset, &offset, U.getDescriptor(),
    VT.localValues, &offset, &offset, VT.getDescriptor(),
    work, &workCount, realWork,
    &info
  );
  for (int i(0); i < A.lens[0]; ++i) {
    LOG(1, "SVG") << "sigma[" << i << "]=" << sigma[i] << std::endl;
  }
  U.write(inverse);
  delete[] sigma;
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

