#include <math/IterativePseudoInverse.hpp>

#include <math/MathFunctions.hpp>
#include <math/RandomTensor.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <complex>
#include <array>
#include <limits>

using namespace cc4s;
using namespace CTF;


template <typename F>
IterativePseudoInverse<F>::IterativePseudoInverse(
  Tensor<F> const &matrix_, F accuracy
):
  matrix(matrix_),
  square(
    2, std::array<int,2>{{matrix_.lens[0], matrix_.lens[0]}}.data(),
    std::array<int,2>{{NS,NS}}.data(), *matrix_.wrld
  ),
  inverse(
    2, std::array<int,2>{{matrix_.lens[0], matrix_.lens[1]}}.data(),
    std::array<int,2>{{NS,NS}}.data(), *matrix_.wrld
  ),
  alpha()
{
  Tensor<F> conjugate(
    2, std::array<int,2>{{matrix.lens[1], matrix.lens[0]}}.data(),
    std::array<int,2>{{NS,NS}}.data(), *matrix.wrld
  );
  Univar_Function<F> fConj(&conj<F>);
  conjugate.sum(1.0,matrix,"ij", 0.0,"ji",fConj);
  square["ij"] = matrix["ik"] * conjugate["kj"];
  Univar_Function<F> fAbs(&abs<F>);
  Vector<F> rowAbsNorms(square.lens[0], *matrix.wrld);
  rowAbsNorms.sum(1.0,square,"ij", 0.0,"i",fAbs);
  std::vector<F> normValues(rowAbsNorms.lens[0]);
  rowAbsNorms.read_all(normValues.data());
  // [K.L. 11.07.2019] complex infinity has undefined behaviour depending
  // on the compiler. Tested using icc-debug and icc config in complex case
  // and it gives (-0,-0). That's why in complex case it works while in real
  // case it doesn't (max=infinity).
  // Anyway, it doesn't make sense to compare if a number is larger than
  // abs(-infinity). 
  // F max(-std::numeric_limits<F>::infinity());
  // A temporary fix: set max default to double type 0. 
  double max(0.);
  for (int i(0); i < square.lens[0]; ++i) {
    if (abs(normValues[i]) > abs(max)) max = abs(normValues[i]);
  }
  alpha = 1.0/max;
  LOG(4, "PseudoInverse") << "alpha=" << alpha << std::endl;
  inverse["ji"] = abs(alpha) * conjugate["ji"];
  iterateQuadratically(accuracy);
//  iterate(accuracy);
}

template <typename F>
void IterativePseudoInverse<F>::iterate(F accuracy) {
  Scalar<F> s;
  Tensor<F> conjugate(
    2, std::array<int,2>{{matrix.lens[1], matrix.lens[0]}}.data(),
    std::array<int,2>{{NS,NS}}.data(), *matrix.wrld
  );
  Univar_Function<F> fConj(&conj<F>);
  conjugate.sum(1.0,matrix,"ij", 0.0,"ji",fConj);
  Tensor<F> sqr(
    2, std::array<int,2>{{matrix.lens[0], matrix.lens[0]}}.data(),
    std::array<int,2>{{NS,NS}}.data(), *matrix.wrld
  );
  F remainder(1.0), minRemainder(std::numeric_limits<F>::infinity());
  int n(0), nMin(0);
  // TODO: use constants for limits
  // TODO: test rectangular matrices with lens[0]>lens[1] & lens[0]>lens[1]
  while (
    abs(remainder) > abs(accuracy*accuracy) &&
    n-nMin < 100 && n < 10000
  ) {
    sqr["ij"] = -1.0 * inverse["ik"] * matrix["kj"];
    sqr["ii"] += 1.0;
    s[""] = sqr["ij"] * sqr["ij"];
    inverse["ij"] += abs(alpha) * sqr["ik"] * conjugate["kj"];
    remainder = abs(s.get_val());
    if (abs(remainder) < abs(minRemainder)) {
      minRemainder = abs(remainder);
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
void IterativePseudoInverse<F>::iterateQuadratically(F accuracy) {
  Scalar<F> s(*matrix.wrld);
  F remainder(1.0), minRemainder(std::numeric_limits<F>::infinity());
  int n(0), nMin(0);
  // TODO: use constants for limits
  while (
    abs(remainder) > abs(accuracy*accuracy) &&
    n-nMin < 2 && n < 10000
  ) {
    square["ij"] = -1.0 * matrix["ik"] * inverse["kj"];
    square["ii"] += 2.0;
    inverse["ij"] = inverse["ik"] * square["kj"];
    square["ii"] += -1.0;
    s[""] = square["ij"] * square["ij"];
    remainder = abs(s.get_val());
    LOG(4, "PseudoInverse") << "remainder=" << remainder << std::endl;
    if (abs(remainder) < abs(minRemainder)) {
      minRemainder = abs(remainder);
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
Tensor<F> &IterativePseudoInverse<F>::get() {
  return inverse;
}

// instantiate
template
class IterativePseudoInverse<cc4s::Float64>;

// not yet there...
//template
//class IterativePseudoInverse<cc4s::Float128>;

template
class IterativePseudoInverse<cc4s::Complex64>;

//template
//class IterativePseudoInverse<cc4s::Complex128>;

template <typename F>
void IterativePseudoInverse<F>::generateHilbertMatrix(Tensor<F> &m) {
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
  Tensor<F> m(
    2, std::array<int,2>{{5,8}}.data(), std::array<int,2>{{NS,NS}}.data(),
    *world
  );
  {
    generateHilbertMatrix(m);
    IterativePseudoInverse pseudoInverse(m);
    Tensor<F> im(pseudoInverse.get());
    im["ij"] = m["ik"] * im["kj"];
    im["ii"] += -1.0;
    Scalar<F> s(*world);
    s[""] = im["ij"] * im["ij"];
    F n(abs(s.get_val()));
    LOG(3) << n << std::endl;
  }
  {
    DefaultRandomEngine random;
    std::normal_distribution<
      typename ComplexTraits<F>::ExtendedType
    > normalDistribution(
      0.0, 1.0
    );
    setRandomTensor(m, normalDistribution, random);
    IterativePseudoInverse pseudoInverse(m);
    Tensor<F> im(pseudoInverse.get());
    im["ij"] = m["ik"] * im["kj"];
    im["ii"] += -1.0;
    Scalar<F> s(*world);
    s[""] = im["ij"] * im["ij"];
    F n(abs(s.get_val()));
    LOG(3) << n << std::endl;
  }
}

// instantiate
template
void IterativePseudoInverse<cc4s::Float64>::test(World *world);
template
void IterativePseudoInverse<cc4s::Complex64>::test(World *world);


template <typename F>
DryIterativePseudoInverse<F>::DryIterativePseudoInverse(
  DryTensor<F> const &matrix_
):
  matrix(matrix_),
  square(
    2, std::array<int,2>{{matrix_.lens[0], matrix_.lens[0]}}.data(),
    std::array<int,2>{{NS,NS}}.data(), SOURCE_LOCATION
  ),
  inverse(
    2, std::array<int,2>{{matrix_.lens[0], matrix_.lens[1]}}.data(),
    std::array<int,2>{{NS,NS}}.data(), SOURCE_LOCATION
  )
{
  DryTensor<F> conjugate(
    2, std::array<int,2>{{matrix_.lens[0], matrix_.lens[1]}}.data(),
    std::array<int,2>{{NS,NS}}.data(), SOURCE_LOCATION
  );
  DryVector<F> rowAbsNorms(square.lens[0]);
}

template <typename F>
DryTensor<F> &DryIterativePseudoInverse<F>::get() {
  return inverse;
}

// instantiate
template
class DryIterativePseudoInverse<cc4s::Float64>;

template
class DryIterativePseudoInverse<cc4s::Complex64>;

