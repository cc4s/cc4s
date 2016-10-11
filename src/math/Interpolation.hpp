#ifndef INTERPOLATION_DEFINED
#define INTERPOLATION_DEFINED

#include <math/MathFunctions.hpp>
#include <iostream>
namespace cc4s {
  template<typename F=double>
  class Interpolation{
  //This class provides, so far, linear interpolation in 1D, 2D and 3D spaces.
  // Further interpolation methods like cubic spline interpolation will be 
  // added.
  public:
    F Linear(F target, F *v){
    // target is a point between 0.--1.0 (normalised)
    // v is a pointer to an array of size 2 containing the two values of the 
    //two points which encompass the target point.
      return (F)((1.-target) * (v[0]) + target * (v[1]));
    }
   
    F Bilinear(F *target, F *v){
    // Bilinear interpolates in 2D space
    // target is the pointer points to a 2D point (x, y), x, y are normalised
    // v is a pointer points to an array of 4 values of the 4
    // points surrounding (x, y), counterclock wise starting from bottom left.
      F v_prime[2] = { Linear(target[0], &(v[0])),
                       Linear(target[0], &(v[2]))};
      return Linear(target[1], v_prime);
    }
    
    F Trilinear(F *target, F *v){
    /*Trilinear interpolates in 3D space.
      target is the pointer points to a 3D point (x, y, z) (normalised).
      v is the pointer points to an array of 8 values of the 8 points which 
      surround (x, y, z).
      The front face of the cubic lies on the plane x = x0 (constant),
      The back face of the cubic lies on  the plane x =x1(constant),
      and x1<x0. 
      Point 0, 1, 2, 3 lie on the front face and counted from the bottom left 
      cornor counterclock wisely;
      Point 4, 5, 6, 7 lie on the back face and counted from the bottom left
      cornor counterclock wisely.
   */
      F v_prime[2] = { Bilinear(&(target[1]), &(v[0])),
                       Bilinear(&(target[1]), &(v[4]))};
      return Linear(target[0], v_prime);
    }
  };
  template <typename F=double>
  inline std::ostream &operator << (std::ostream &stream, cc4s::Interpolation<F> const &v) {
    stream << v << std::endl;
    return stream;
  }
}
#endif
