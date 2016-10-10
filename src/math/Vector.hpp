#ifndef VECTOR_DEFINED
#define VECTOR_DEFINED

#include <math/MathFunctions.hpp>
#include <iostream>

namespace cc4s {
  template <typename F=double, int D=3>
  class Vector {
  public:
    Vector() {
      for (int d(0); d<D; ++d) {
        coordinate[d] = 0;
      }
    }
    Vector(Vector<F,D> const &v) {
      for (int d(0); d<D; ++d) {
        coordinate[d] = v.coordinate[d];
      }
    }

    Vector<F,D> &operator += (Vector<F,D> const &v) {
      for (int d(0); d<D; ++d) { 
        coordinate[d] += v.coordinate[d];
      }
      return *this;
    }

    Vector<F,D> &operator -= (Vector<F,D> const &v) {
      for (int d(0); d<D; ++d) { 
        coordinate[d] -= v.coordinate[d];
      }
      return *this;
    }

    F dot(Vector<F,D> const &v) {
      F sum(0);
      for (int d(0); d<D; ++d) { 
        sum += coordinate[d] * cc4s::conj(v.coordinate[d]);
      }
      return sum;
    }

    int approximately(Vector<F,D> const &v, const double epsilon = 1e-15) {
      Vector<F,D> u(*this);
      u -= v;
      return std::real(u.dot(u)) < epsilon;
    }

    F coordinate[D];
  };

  template <typename F=double, int D=3>
  inline std::ostream &operator << (std::ostream &stream, cc4s::Vector<F,D> const &v) {
    for (int d(0); d<D-1; ++d) {
      stream << v.coordinate[d] << ",";
    }
    stream << v.coordinate[D-1];
    return stream;
  }
}

#endif

