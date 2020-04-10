#ifndef INTEGRATION_DEFINED
#define INTEGRATION_DEFINED

namespace cc4s {
  template<typename F, typename Real, typename Method>
  Real integrate(F f, Real a, Real b, size_t steps, Method m) {
    Real s(0);
    Real h((b-a)/steps);
    for (size_t i(0); i <= steps; ++i) {
      s += m(f, (a*(steps-i) + b*i)/steps, h);
    }
    return h*s;
  }

  class Trapezium {
  public:
    template<typename F, typename Real>
    Real operator()(F f, Real x, Real h) const {
      return (f(x) + f(x+h))/2;
    }
  };
 
  class Simpson {
  public:
    template<typename F, typename Real>
    Real operator()(F f, Real x, Real h) const {
      return (f(x) + 4*f(x+h/2) + f(x+h))/6;
    }
  };
}
#endif
