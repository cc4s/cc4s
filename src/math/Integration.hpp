#ifndef INTEGRATION_DEFINED
#define INTEGRATION_DEFINED

#include <math/MathFunctions.hpp>
#include <iostream>
#include <string>
namespace cc4s {
  template<typename F, typename Real, typename Method>
  Real integrate(F f, Real a, Real b, int steps, Method m){
    double s = 0;
    double h = (b-a)/steps;
    for (int i = 0; i < steps; ++i)
      s += m(f, a + h*i, h);
    return h*s;
  }

  class trapezium{
  public:
    template<typename F, typename Real>
    double operator()(F f, Real x, Real h) const{
      return (f(x) + f(x+h))/2;
    }
  };
 
  class simpson{
  public:
    template<typename F, typename Real>
    double operator()(F f, Real x, Real h) const{
      return (f(x) + 4*f(x+h/2) + f(x+h))/6;
    }
  };
}
#endif
