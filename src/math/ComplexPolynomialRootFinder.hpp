/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COMPLEX_POLYNOMIAL_ROOT_FINDER_DEFINED
#define COMPLEX_POLYNOMIAL_ROOT_FINDER_DEFINED

#include <math/Complex.hpp>

namespace util {
  template <typename T>
  class Polynomial {
  public:
    Polynomial(
      Polynomial const &p
    ): a(new T[p.n+1]), n(p.n) {
      for (int k(0); k < n; ++k) {
        a[k] = p.a[k];
      }
    }

    Polynomial(T const *a_, int const n_): a(new T[n_ + 1]), n(n_) {
      for (int k(0); k <= n; ++k) {
        a[k] = a_[k];
      }
    }
    ~Polynomial() {
      delete []a;
      a = nullptr;
    }

    int degree() const {
      return n;
    }

    T at(T const x) const {
      T p(a[n]);
      for (int k(n-1); k >= 0; --k) {
        p = p*x + a[k];
      }
      return p;
    }

    T derivativeAt(T const x) const {
      T p(n*a[n]);
      for (int k(n-1); k >= 1; --k) {
        p = p*x + T(k)*a[k];
      }
      return p;
    }

  protected:
    T *a;
    int n;
  };

  template <typename T>
  class ScaledPolynomial {
  public:
    ScaledPolynomial(
      ScaledPolynomial const &p
    ): a(new T[p.n]), n(p.n), an(p.an), shift(p.shift) {
      for (int k(0); k < n; ++k) {
        a[k] = p.a[k];
      }
    }

    ScaledPolynomial(
      T const *a_, int const n_
    ): a(new T[n_]), n(n_), an(a_[n_]), shift(0) {
      // scale the polynomial in y
      for (int k(0); k < n; ++k) {
        a[k] = a_[k] / an;
        int sk;
        std::frexp(std::abs(a[k]), &sk);
        shift = std::max(shift, sk/(n-k));
      }
      // scale the polynomial in x
      for (int k(0); k < n; ++k) {
        a[k] *= std::ldexp(1.0, shift*(k-n));;
      }
    }
    ~ScaledPolynomial() {
      delete []a;
      a = nullptr;
    }

    int degree() const {
      return n;
    }

    T at(T const x) const {
      T p(1);
      for (int k(n-1); k >= 0; --k) {
        p = p*x + a[k];
      }
      return p;
    }

    T derivativeAt(T const x) const {
      T p(n);
      for (int k(n-1); k >= 1; --k) {
        p = p*x + T(k)*a[k];
      }
      return p;
    }

    double getScale() const {
      return std::ldexp(1.0, shift);
    }

  protected:
    T *a;
    int n;
    T an;
    int shift;
  };

  /**
  * \brief Finds complex roots of a polynomial of arbitrary ordrer.
  */
  class ComplexPolynomialRootFinder {
  public:
    ComplexPolynomialRootFinder(
      ScaledPolynomial<Complex<>> const * const p_
    ): p(p_), roots(new Complex<>[p_->degree()]) {
      int rnd(1);
      // initialize the roots randomly wihtin the unit square [-1,1]^2
      // this is an upper bound for the roots of a scaled polynomial
      for (int k(0); k < p->degree(); ++k) {
        rnd = (rnd*21839 + 3643) % 32749;
        double re(2.0 * rnd/32749 - 1.0);
        rnd = (rnd*21839 + 3643) % 32749;
        double im(2.0 * rnd/32749 - 1.0);
        roots[k] = Complex<>(re,im);
      }
    }

    void findRoots(double const accuracy = 1e-10) {
      Complex<> steps[p->degree()];
      double epsilon;
      do {
        for (int k(0); k < p->degree(); ++k) {
          Complex<> r(
            p->at(roots[k]) / p->derivativeAt(roots[k])
          );
          Complex<> repulsion(0);
          for (int j(0); j < p->degree(); ++j) {
            if (j != k) repulsion +=
              Complex<>(1) / (roots[k] - roots[j]);
          }
          steps[k] = r / (Complex<>(1) - r * repulsion);
        }
        epsilon = 0.0;
        for (int k(0); k < p->degree(); ++k) {
          epsilon += steps[k].real() * steps[k].real();
          epsilon += steps[k].imag() * steps[k].imag();
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
    ScaledPolynomial<Complex<>> const *p;
    Complex<> *roots;
  };
}

#endif

