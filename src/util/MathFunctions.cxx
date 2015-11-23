#include "MathFunctions.hpp"
#include <cmath>

// univariate functions
double MathFunctions::sqrt(double const x) {
  return std::sqrt(x);
}


template <typename F>
F MathFunctions::conj(F const z) {
  return std::conj(z);
}
// instantiate:
template double MathFunctions::conj(double const x);
template complex MathFunctions::conj(complex const x);

template <typename F>
F MathFunctions::abs(F const z) {
  return std::abs(z);
}
// instantiate:
template double MathFunctions::abs(double const x);
template complex MathFunctions::abs(complex const x);

// bivariate functions
double MathFunctions::divide(double const x, double const y) {
  return x / y;
}

