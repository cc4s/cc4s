#ifndef MATH_FUNCTIONS_DEFINED
#define MATH_FUNCTIONS_DEFINED

#include <util/Complex.hpp>

/**
 * \deprecated std::abs<double>, ... can be instantiated for the same purpose
 */ 
class MathFunctions {
public:
  // univariate functions
  static double sqrt(double const x);
  template <typename F>
  static F conj(F const z);
  template <typename F>
  static F abs(F const z);

  // bivariate functions
  static double divide(double const x, double const y);
};

#endif
