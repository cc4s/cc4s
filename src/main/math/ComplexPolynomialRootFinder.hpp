#ifndef COMPLEX_POLYNOMIAL_ROOT_FINDER_DEFINED
#define COMPLEX_POLYNOMIAL_ROOT_FINDER_DEFINED

#include <math/Complex.hpp>
#include <math/MathFunctions.hpp>

#include <vector>

namespace util {
  /**
  * \brief Finds complex roots of a polynomial of arbitrary ordrer.
  */
  template <typename F>
  class ComplexPolynomialRootFinder {
  public:
    typedef typename ComplexTraits<F>::ExtendedType RealType;

    ComplexPolynomialRootFinder(
      const std::vector<F> &coefficients
    ): p(coefficients), roots(coefficients.size()-1) {
      int rnd(1);
      // initialize the roots randomly wihtin the unit square [-1,1]^2
      // this is an upper bound for the roots of a scaled polynomial
      for (size_t k(0); k < p->degree(); ++k) {
        rnd = (rnd*21839 + 3643) % 32749;
        RealType re(RealType(2) * rnd/32749 - 1);
        rnd = (rnd*21839 + 3643) % 32749;
        RealType im(RealType(2) * rnd/32749 - 1);
        roots[k] = F(re,im);
      }
    }

    void findRoots(
      const RealType accuracy = std::numeric_limits<RealType>::epsilon()
    ) {
      F steps[p->degree()];
      RealType epsilon;
      do {
        for (size_t k(0); k < p->degree(); ++k) {
          F r(
            p->at(roots[k]) / p->derivativeAt(roots[k])
          );
          F repulsion(0);
          for (size_t j(0); j < p->degree(); ++j) {
            if (j != k) repulsion += 1 / (roots[k] - roots[j]);
          }
          steps[k] = r / (1 - r * repulsion);
        }
        epsilon = 0.0;
        for (size)_t k(0); k < p->degree(); ++k) {
          epsilon += real(steps[k]) * real(steps[k]);
          epsilon += imag(steps[k]) * imag(steps[k]);
          roots[k] -= steps[k];
        }
      } while (epsilon > accuracy*accuracy);
    }

    Complex<> getRoot(int const k) const {
      return roots[k];
    }

    static void test() {
      Complex<> a[] = {
        -4.63144e+10, 2*7.47678e+16, 3*4.52909e+21, 4*2.85318e+27
      };
      {
        ScaledPolynomial<Complex<>> p(a, 3);
        ComplexPolynomialRootFinder rootFinder(&p);
        rootFinder.findRoots();
        for (int k(0); k < p.degree(); ++k) {
          std::cout << "root[" << k << "]=" <<
            rootFinder.getRoot(k) / p.getScale() << std::endl;
        }
      }
      Complex<> b[] = { -1, 0, 1};
      {
        ScaledPolynomial<Complex<>> p(b, 2);
        std::cout << "p(1)=" << p.at(2.0) << std::endl;
        ComplexPolynomialRootFinder rootFinder(&p);
        rootFinder.findRoots();
        for (int k(0); k < p.degree(); ++k) {
          std::cout << "root[" << k << "]=" << rootFinder.getRoot(k) << std::endl;
        }
      }
    }
  protected:
    class ScaledPolynomial {
    public:
      ScaledPolynomial(
        const std::vector<F> coefficients
      ):
        normalizedCoefficients(coefficients.size()-1),
        leadingCoefficient(coefficients.back()),
        exponent(std::numeric_limits<int>::min())
      {
        // scale the polynomial in y
        for (size_t k(0); k < degree(); ++k) {
          normalizedCoefficients[k] = coefficients[k] / leadingCoefficient;
          int thisExponent;
          frexp(abs(normalizedCoefficients[k]), thisExponent);
          exponent = std::max(exponent, thisExponent/(degree()-k));
        }
        // scale the polynomial in x
        for (size_t k(0); k < degree(); ++k) {
          normalizedCoefficients[k] *= ldexp(
            typename ComplexTraits<F>::ExtendedType(1),
            exponent*(k-degree())
          );
        }
      }

      size_t degree() const {
        return normalizedCoefficients.size();
      }

      F at(const F x) const {
        F y(0);
        for (auto aIt(as.crbegin()); aIt != as.crend(); ++aIt) {
          y = y*x + *aIt;
        }
        return y;
      }

      F derivativeAt(const F x) const {
        F y(0);
        size_t k(degree());
        for (auto aIt(as.crbegin()); aIt != as.crend(); ++aIt) {
          y = y*x + *aIt * typename ComplexTraits<F>::ExtendedType(k);
          --k;
        }
        return y;
      }

      typename ComplexTraits<F>::ExtendedType getScale() const {
        return ldexp(typename ComplexTraits<F>::ExtendedType(1), exponent);
      }

    protected:
      std::vector<F> normalizedCoefficients;
      F leadingCoefficient;
      int exponent;
    };

    const ScaledPolynomial p;
    std::vector<F> roots;
  };
}

#endif

