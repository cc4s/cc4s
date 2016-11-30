#include <math/IterativePseudoInverse.hpp>

#include <math/MathFunctions.hpp>
#include <math/RandomTensor.hpp>
#include <util/Log.hpp>
#include <complex>
#include <limits>

using namespace cc4s;
using namespace CTF;


template <typename F>
IterativePseudoInverse<F>::IterativePseudoInverse(
  Matrix<F> const &matrix_, double accuracy
):
  matrix(matrix_),
  square(matrix_.lens[0], matrix_.lens[0], NS, *matrix_.wrld),
  inverse(matrix_.lens[1], matrix_.lens[0], NS, *matrix_.wrld),
  alpha()
{
  Matrix<F> conjugate(matrix.lens[1], matrix.lens[0], NS, *matrix.wrld);
  Univar_Function<F> fConj(&conj<F>);
  conjugate.sum(1.0,matrix,"ij", 0.0,"ji",fConj);
  Matrix<F> square(matrix.lens[0], matrix.lens[0], NS, *matrix.wrld);
  square["ij"] = matrix["ik"] * conjugate["kj"];
  Univar_Function<F> fAbs(&abs<F>);
  Vector<F> rowAbsNorms(square.lens[0], *matrix.wrld);
  rowAbsNorms.sum(1.0,square,"ij", 0.0,"i",fAbs);
  F *normValues(new F[rowAbsNorms.lens[0]]);
  rowAbsNorms.read_all(normValues);
  double max(-std::numeric_limits<double>::infinity());
  for (int i(0); i < square.lens[0]; ++i) {
    if (std::real(normValues[i]) > max) max = std::real(normValues[i]);
  }
  alpha = 1.0/max;
  LOG(4, "PseudoInverse") << "alpha=" << alpha << std::endl;
  inverse["ji"] = alpha * conjugate["ji"];
  iterateQuadratically(accuracy);
//  iterate(accuracy);
}

template <typename F>
void IterativePseudoInverse<F>::iterate(double accuracy) {
  Scalar<F> s;
  Matrix<F> conjugate(matrix.lens[1], matrix.lens[0], NS, *matrix.wrld);
  Univar_Function<F> fConj(&conj<F>);
  conjugate.sum(1.0,matrix,"ij", 0.0,"ji",fConj);
  Matrix<F> sqr(matrix.lens[1], matrix.lens[1], *matrix.wrld);  double remainder(1.0), minRemainder(std::numeric_limits<double>::infinity());
  int n(0), nMin(0);
  // TODO: use constants for limits
  while (remainder > accuracy*accuracy && n-nMin < 100 && n < 10000) {

    sqr["ij"] = -1.0 * inverse["ik"] * matrix["kj"];
    sqr["ii"] += 1.0;
    Bivar_Function<F> fRealDot(&realDot<F>);
    s.contract(1.0, sqr,"ij", sqr,"ij", 0.0,"", fRealDot);
    inverse["ij"] += alpha * sqr["ik"] * conjugate["kj"];
    remainder = std::real(s.get_val());
    if (remainder < minRemainder) {
      minRemainder = remainder;
      nMin = n;
    }
    ++n;
  }
  if (n >= 10000) {
    // failed to convege
    LOG(0, "PseudoInverse") << "failed to converge, remainder=" << remainder
      << ", minRemainder=" << minRemainder << std::endl;
    throw new EXCEPTION("Failed to converge iterative pseudo inverse.");
  }
}

template <typename F>
void IterativePseudoInverse<F>::iterateQuadratically(double accuracy) {
  Scalar<F> s(*matrix.wrld);
  double remainder(1.0), minRemainder(std::numeric_limits<double>::infinity());
  int n(0), nMin(0);
  // TODO: use constants for limits
  while (remainder > accuracy*accuracy && n-nMin < 2 && n < 10000) {
    square["ij"] = -1.0 * matrix["ik"] * inverse["kj"];
    square["ii"] += 2.0;
    inverse["ij"] = inverse["ik"] * square["kj"];
    square["ii"] += -1.0;
    Bivar_Function<F> fRealDot(&realDot<F>);
    s.contract(1.0, square,"ij", square,"ij", 0.0,"", fRealDot);
    remainder = std::real(s.get_val());
    LOG(4, "PseudoInverse") << "remainder=" << remainder << std::endl;
    if (remainder < minRemainder) {
      minRemainder = remainder;
      nMin = n;
    }
    ++n;
  }
  if (n >= 10000) {
    // failed to convege
    LOG(0, "PseudoInverse") << "failed to converge, remainder=" << remainder
      << ", minRemainder=" << minRemainder << std::endl;
    throw new EXCEPTION("Failed to converge iterative pseudo inverse.");
  }
}

template <typename F>
Matrix<F> &IterativePseudoInverse<F>::get() {
  return inverse;
}

// instantiate
template
IterativePseudoInverse<double>::IterativePseudoInverse(
  Matrix<double> const &matrix, double accuracy
);
template
Matrix<double> &IterativePseudoInverse<double>::get();

template
IterativePseudoInverse<complex>::IterativePseudoInverse(
  Matrix<complex> const &matrix, double accuracy
);
template
Matrix<complex> &IterativePseudoInverse<complex>::get();


template <typename F>
void IterativePseudoInverse<F>::generateHilbertMatrix(Matrix<F> &m) {
  int64_t indicesCount, *indices;
  F *values;
  m.read_local(&indicesCount, &indices, &values);
  for (int64_t l(0); l < indicesCount; ++l) {
    int i = int(l % m.lens[0]);
    int j = int(l / m.lens[0]);
    values[l] = 1.0 / (i+j+1);
  }
  m.write(indicesCount, indices, values);
  free(indices); free(values);
}

template <typename F>
void IterativePseudoInverse<F>::test(World *world) {
  Matrix<F> m(5, 8, NS, *world);
  {
    generateHilbertMatrix(m);
    IterativePseudoInverse pseudoInverse(m);
    Matrix<F> im(pseudoInverse.get());
    im["ij"] = m["ik"] * im["kj"];
    im["ii"] += -1.0;
    Scalar<F> s(*world);
    s[""] = im["ij"] * im["ij"];
    double n(std::real(s.get_val()));
    LOG(3) << n << std::endl;
  }
  {
    setRandomTensor(m);
    IterativePseudoInverse pseudoInverse(m);
    Matrix<F> im(pseudoInverse.get());
    im["ij"] = m["ik"] * im["kj"];
    im["ii"] += -1.0;
    Scalar<F> s(*world);
    s[""] = im["ij"] * im["ij"];
    double n(std::real(s.get_val()));
    LOG(3) << n << std::endl;
  }
}

// instantiate
template
void IterativePseudoInverse<double>::test(World *world);
template
void IterativePseudoInverse<complex>::test(World *world);


template <typename F>
DryIterativePseudoInverse<F>::DryIterativePseudoInverse(
  DryMatrix<F> const &matrix_
):
  matrix(matrix_),
  square(matrix_.lens[0], matrix_.lens[0], NS),
  inverse(matrix_.lens[1], matrix_.lens[0], NS)
{
  DryMatrix<F> conjugate(matrix.lens[1], matrix.lens[0], NS);
  DryMatrix<F> square(matrix.lens[0], matrix.lens[0], NS);
  DryVector<F> rowAbsNorms(square.lens[0]);
}

template <typename F>
DryMatrix<F> &DryIterativePseudoInverse<F>::get() {
  return inverse;
}

// instantiate
template
DryIterativePseudoInverse<double>::DryIterativePseudoInverse(
  DryMatrix<double> const &matrix
);
template
DryMatrix<double> &DryIterativePseudoInverse<double>::get();

template
DryIterativePseudoInverse<complex>::DryIterativePseudoInverse(
  DryMatrix<complex> const &matrix
);
template
DryMatrix<complex> &DryIterativePseudoInverse<complex>::get();

