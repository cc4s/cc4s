#include "MathFunctions.hpp"
#include <cmath>

// univariate functions
double MathFunctions::sqrt(double const x) {
  return std::sqrt(x);
}

complex MathFunctions::conj(complex const z) {
  return std::conj(z);
}

// bivariate functions
double MathFunctions::divide(double const x, double const y) {
  return x / y;
}

