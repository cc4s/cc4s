#include <util/IterativePseudoInverter.hpp>

#include <util/MathFunctions.hpp>
#include <util/Log.hpp>
#include <limits>

using namespace CTF;

// TODO: place in proper file
template <typename F>
void dumpMatrix(Matrix<F> &m) {
  F *values(new F[m.lens[0]*m.lens[1]]);
  m.read_all(values);
  for (int i(0); i < m.lens[0]; ++i) {
    for (int j(0); j < m.lens[1]; ++j) {
      LOG(3) << values[i+j*m.lens[0]];
    }
    LOG(3) << std::endl;
  }
}

template <typename F>
IterativePseudoInverter<F>::IterativePseudoInverter(
  Matrix<F> const &matrix_
):
  matrix(matrix_),
  inverse(matrix_.lens[1], matrix_.lens[0], *matrix.wrld),
  square(matrix_.lens[0], matrix_.lens[0], *matrix.wrld)
{
  Matrix<F> conjugate(matrix.lens[1], matrix.lens[0], *matrix.wrld);
  Univar_Function<F> fconj(&MathFunctions::conj<F>);
  conjugate.sum(1.0,matrix,"ij", 0.0,"ji",fconj);
  square["ij"] = matrix["ik"] * conjugate["kj"];
//  Univar_Function<F> fabs(&MathFunctions::abs<F>);
  Univar_Function<F> fabs(&std::abs<F>);
  Vector<F> rowAbsNorms(square.lens[0], *matrix.wrld);
  rowAbsNorms.sum(1.0,square,"ij", 0.0,"i",fabs);
  F *normValues(new F[square.lens[0]]);
  rowAbsNorms.read_all(normValues);
  double max(-std::numeric_limits<double>::infinity());
  for (int i(0); i < square.lens[0]; ++i) {
    if (std::real(normValues[i]) > max) max = std::real(normValues[i]);
  }
  inverse["ji"] = (2.0/max) * conjugate["ji"];
}

template <typename F>
void IterativePseudoInverter<F>::invert(
  Matrix<F> &pseudoInverse, double accuracy
) {
  Scalar<F> s;
  for (int n(0); n < 10000; ++n) {
    square["ij"] = -1.0 * matrix["ik"] * inverse["kj"];
    square["ii"] += 2.0;
    inverse["ij"] = inverse["ik"] * square["kj"];
    square["ii"] += -1.0;
    s[""] = square["ij"] * square["ij"];
    if (std::real(s.get_val()) < accuracy*accuracy) {
      LOG(4) << "  converged after " << n << " steps" << std::endl;
      pseudoInverse["ij"] = inverse["ij"];
      return;
    }
  }
  double rest(std::real(s.get_val()));
  // failed to convege
  LOG(4) << " failed to converge, rest=" << rest << std::endl;
  dumpMatrix(inverse);
}

// instantiate
template
IterativePseudoInverter<double>::IterativePseudoInverter(
  Matrix<double> const &matrix
);
template
void IterativePseudoInverter<double>::invert(Matrix<double> &pseudoInverse, double accuracy);

template
IterativePseudoInverter<complex>::IterativePseudoInverter(
  Matrix<complex> const &matrix
);
template
void IterativePseudoInverter<complex>::invert(Matrix<complex> &pseudoInverse, double accuracy);



template <typename F>
void IterativePseudoInverter<F>::test(World *world) {
  // generate Hilbert Matrix
  Matrix<F> m(10, 10, NS, *world);
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
  IterativePseudoInverter pseudoInverter(m);
  Matrix<F> im(m.lens[1], m.lens[0], NS, *world);
  pseudoInverter.invert(im);
  dumpMatrix(im);
  im["ij"] = m["ik"] * im["kj"];
  im["ii"] += -1.0;
  Scalar<F> s(*world);
  s[""] = im["ij"] * im["ij"];
  double n(std::real(s.get_val()));
  LOG(3) << n << std::endl;
}

// instantiate
template
void IterativePseudoInverter<double>::test(World *world);
template
void IterativePseudoInverter<complex>::test(World *world);
