#ifndef SAMPLED_VARIABLE_DEFINED
#define SAMPLED_VARIABLE_DEFINED

#include <math/Complex.hpp>

namespace cc4s {
  template <typename F=double>
  class SampledVariable {
  protected:
  public:
    SampledVariable(): s(0), s2(0), n(0) {
    }

    void addSample(F x) {
      s += x;
      s2 += absSqr(x);
      ++n;
    }

    F getMean() {
      return s / static_cast<double>(n);
    }

    F getVariance() {
      F m(getMean());
      return ( s2 - n*absSqr(m) ) / (n-1);
    }

    F getMeanVariance() {
      return getVariance() / static_cast<double>(n);
    }

    F getMeanStdDeviation() {
      return std::sqrt(getMeanVariance());
    }

    F s;
    double s2;
    uint64_t n;
  };
}

#endif

