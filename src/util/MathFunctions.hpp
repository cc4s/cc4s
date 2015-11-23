#ifndef MATH_FUNCTIONS_DEFINED
#define MATH_FUNCTIONS_DEFINED

#include <Complex.hpp>

class MathFunctions {
public:
  // univariate functions
  static double sqrt(double const x);
  static complex conj(complex const z);

  // bivariate functions
  static double divide(double const x, double const y);
};

#endif

