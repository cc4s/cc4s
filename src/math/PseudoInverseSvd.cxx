#include <math/PseudoInverseSvd.hpp>

#include <math/MathFunctions.hpp>
#include <util/BlacsWorld.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <extern/ScaLapack.hpp>
#include <util/Log.hpp>

using namespace cc4s;
using namespace CTF;


template <typename F>
PseudoInverseSvd<F>::PseudoInverseSvd(
  Matrix<F> const &matrix_
):
  A(matrix_)
{
}

template <typename F>
Matrix<F> &PseudoInverseSvd<F>::get() {
  BlacsWorld blacsWorld(A.wrld->rank, A.wrld->np);
  ScaLapackMatrix<F> ScaA(A, &blacsWorld);
  ScaLapackMatrix<F> ScaU(ScaA);
  ScaLapackMatrix<F> ScaVT(ScaA);
  int offset(1);
  double *sigma(new double[A.lens[0]]);
  F optimalWork;
  double optimalRealWork;
  int info;
  int workCount(-1);
  pgesvd(
    "V", "V",
    &ScaA.lens[0], &ScaA.lens[1],
    ScaA.localValues, &offset, &offset, ScaA.getDescriptor(),
    sigma,
    ScaU.localValues, &offset, &offset, ScaU.getDescriptor(),
    ScaVT.localValues, &offset, &offset, ScaVT.getDescriptor(),
    &optimalWork, &workCount, &optimalRealWork,
    &info
  );
  LOG(4, "SVD") << "work size=" << optimalWork << "," << optimalRealWork << std::endl;
  workCount = static_cast<int>(std::real(optimalWork)+0.5);
  F *work(new F[workCount]);
  double *realWork(new double[static_cast<int64_t>(optimalRealWork+0.5)]);
  pgesvd(
    "V", "V",
    &ScaA.lens[0], &A.lens[1],
    ScaA.localValues, &offset, &offset, ScaA.getDescriptor(),
    sigma,
    ScaU.localValues, &offset, &offset, ScaU.getDescriptor(),
    ScaVT.localValues, &offset, &offset, ScaVT.getDescriptor(),
    work, &workCount, realWork,
    &info
  );
  delete[] work; delete[] realWork;
  Vector<F> S(A.lens[0], *A.wrld, "Sigma");
  int localSigmaCount(A.wrld->rank == 0 ? ScaA.lens[0] : 0);
  F *sigmaValues(new F[localSigmaCount]);
  int64_t *sigmaIndices(new int64_t[localSigmaCount]);
  // TODO: invert singular values
  for (int64_t i(0); i < localSigmaCount; ++i) {
    sigmaIndices[i] = i;
    sigmaValues[i] = sigma[i];
  }
  S.write(localSigmaCount, sigmaIndices, sigmaValues);
  Matrix<F> U(A);
  Matrix<F> VT(A);
  ScaU.write(U);
  ScaVT.write(VT);
  delete[] sigmaIndices; delete[] sigmaValues;
  delete[] sigma;

  U["ij"] = U["ij"] * S["j"];
  A["ij"] = U["ik"] * VT["kj"];
  return A;
}

// instantiate
template
PseudoInverseSvd<double>::PseudoInverseSvd(
  Matrix<double> const &matrix
);
template
Matrix<double> &PseudoInverseSvd<double>::get();

template
PseudoInverseSvd<complex>::PseudoInverseSvd(
  Matrix<complex> const &matrix
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

