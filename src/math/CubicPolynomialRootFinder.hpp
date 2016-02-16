/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CUBIC_POLYNOMIAL_ROOT_FINDER_DEFINED
#define CUBIC_POLYNOMIAL_ROOT_FINDER_DEFINED

// FIXME: write test cases

namespace util {
  // TODO: find approprate place
  template <typename T>
  int sgn(T value) {
    return (T(0) < value) - (value < T(0));
  }

  /**
  * \brief finds minimum and maximum of a polynomial of 4th order.
  * The algorithm follows the suggestions of Deiters and Macias-Salinas
  * http://pubs.acs.org/doi/pdf/10.1021/ie4038664.
  */
  class CubicPolynomialRootFinder {
  public:
    CubicPolynomialRootFinder(
      double const a0_, double const a1_, double const a2_, double const a3_,
      double const epsilon_ = 1e-10
    ): a0(a0_), a1(a1_), a2(a2_), a3(a3_), epsilon(epsilon_) {
    }
    /**
    * \brief find the roots of the polynomial in any order.
    * NOTE: the case that a3 may be zero is not treated.
    * \returns the number of real roots found.
    */
    int findRoots(double roots[3]) {
      // normalize
      a0 /= a3; a1 /= a3; a2 /= a3;
      // scale to make the largest coefficient approximately +/-1
      int s0, s1, s2;
      frexp(a0, &s0); frexp(a1, &s1); frexp(a2, &s2);
      shift = std::max(s0/3, std::max(s1<<1, s2));
      a0 = ldexp(a0,-3*shift); a1 = ldexp(a1,-2*shift); a2 = ldexp(a2,-shift);
      // NOTE that all roots must be scaled by 2^shift to arrive at the
      // roots of the original polynomial
      int rootsCount = findScaledRoots(roots);
      // unscale the roots
      for (int i(0); i < rootsCount; ++i) {
        roots[i] = ldexp(roots[i], shift);
      }
      // unscale the coefficients for later evaluateAt calls
      a0 = ldexp(a0,+3*shift); a1 = ldexp(a1,+2*shift); a2 = ldexp(a2,+shift);
      return rootsCount;
    }

    double evaluateAt(double const x) {
      double r(x+a2);
      r = r*x+a1;
      return r*x + a0;
    }

    double evaluateDerivativeAt(double const x) {
      double r(3.0*x+2.0*a2);
      return r*x+a1;
    }

    double evaluateDerivative2At(double const x) {
      return 6.0*x+a2;
    }

    static void test() {
//      double a0(392317.0);
      double a1(-4.63144e+10);
      double a2(7.47678e+16);
      double a3(4.52909e+21);
      double a4(2.85318e+27);
      CubicPolynomialRootFinder rootFinder(
        a1, 2.0*a2, 3.0*a3, 4.0*a4
      );
      double roots[3];
      int rootsCount(rootFinder.findRoots(roots));
      for (int i(0); i < rootsCount; ++i) {
        std::cout << "r=" << roots[i] << std::endl;
      }
    }
  protected:
    double a0, a1, a2, a3;
    double epsilon;
    int shift;

    int findScaledRoots(double roots[3]) {
      if (std::abs(a0) < epsilon) {
        // one root is zero
        roots[0] = 0.0;
        int n = findQuadraticRoots(a1, a2, &roots[1]);
        return n+1;
      }
      // calculate the inflection point
      double xInf(-a2/3.0), yInf(evaluateAt(xInf));

      if (std::abs(yInf) < epsilon) {
        // the inflection point is a root
        roots[0] = xInf;
        double b0, b1;
        deflate(xInf, b0, b1);
        int n = findQuadraticRoots(b0, b1, &roots[1]);
        return n+1;
      }
      double d(a2*a2 - 3*a1);
      double root;
      if (d <= -epsilon) {
        // initialize the iteration with xInf
        root = xInf;
      } else if (d < epsilon) {
        // the slope is zero
        roots[0] = xInf - std::cbrt(yInf);
        return 1;
      } else {
        // otherwise, initialize as follows
        root = xInf - 2.0*sgn(yInf)/3.0*std::sqrt(d);
      }

      double dx;
      do {
        double y(evaluateAt(root));
        double dy(evaluateDerivativeAt(root));
        double ddy(evaluateDerivative2At(root));
        dx = y*dy/(dy*dy-0.5*y*ddy);
        root -= dx;
      } while (std::abs(dx) > epsilon*std::abs(root));
      roots[0] = root;

      if (d <= epsilon) {
        return 1;
      } else {
        double b0, b1;
        deflate(root, b0, b1);
        int n = findQuadraticRoots(b0, b1, &roots[1]);
        return n + 1;
      }
    }

    int findQuadraticRoots(double const b0, double const b1, double roots[2]) {
      double d(0.25*b1*b1 - b0);
      if (d <= -epsilon) {
        return 0;
      } else if (d < epsilon) {
        roots[0] = roots[1] = -0.5*b1;
      } else  {
        d = sqrt(d);
        roots[0] = -0.5*b1 - d;
        roots[1] = -0.5*b1 + d;
      }
      return 2;
    }

    void deflate(double const root, double &b0, double &b1) {
      b1 = a2 + root;
      b0 = b1*root + a1;
    }
  };
};

#endif

