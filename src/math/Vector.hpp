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
    
    Vector<F,D> &operator / (F r){
      for (int d(0); d<D; ++d){
        coordinate[d] = coordinate[d] / r;
      }
      return *this;
    }

    Vector<F,D> &operator * (F r){
      for (int d(0); d<D; ++d){
        coordinate[d] = coordinate[d] * r;
      }
      return *this;
    }

    // Specialization for D=3
    Vector<F,3> cross(Vector<F,3> const &v) {
      Vector<F,3> t;
      t.coordinate[0] = coordinate[1]*v.coordinate[2] - 
                          coordinate[2]*v.coordinate[1];
      t.coordinate[1] = coordinate[2]*v.coordinate[0] - 
                          coordinate[0]*v.coordinate[2];
      t.coordinate[2] = coordinate[0]*v.coordinate[1] - 
                          coordinate[1]*v.coordinate[0];
      return t;
    }

    F dot(Vector<F,D> const &v) {
      F sum(0);
      for (int d(0); d<D; ++d) { 
        sum += coordinate[d] * cc4s::conj(v.coordinate[d]);
      }
      return sum;
    }

    int approximately(Vector<F,D> const &v, const double epsilon = 1e-10) const {
      Vector<F,D> u(*this);
      u -= v;
      return std::real(u.dot(u)) < epsilon;
    }

    double distance(Vector<F, D> const &v) const {
      Vector<F,D> u(*this);
      u -= v;
      return std::real(u.dot(u)); 
    }

    double length() const {
      Vector<F,D> u(*this);
      return std::sqrt(std::real(u.dot(u))); 
    }

    F operator [](int d) const {
      return coordinate[d];
    }

    F at(int d) const {
      return coordinate[d];
    }

    F &at(int d) {
      return coordinate[d];
    }

    F &operator [](int d) {
      return coordinate[d];
    }
    
    static bool sortByLength(Vector<F, D> const &n, Vector<F,D> const &m){
      return n.length() < m.length();
    }

    bool operator < (Vector<F, D>  const &v) const {
      const double epsilon(1e-10);
      for (int d(0); d<D; ++d) {
        if (coordinate[d] < v.coordinate[d] - epsilon) {
          return true;
        } else if (coordinate[d] > v.coordinate[d] + epsilon) {
          return false;
        }
      }
      return false;
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
  };

}

#endif

