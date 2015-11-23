#ifndef ITERATIVE_PSEUDO_INVERTER_DEFINED
#define ITERATIVE_PSEUDO_INVERTER_DEFINED

#include <util/Complex.hpp>
#include <ctf.hpp>

template <typename F>
class IterativePseudoInverter {
public:
  IterativePseudoInverter(CTF::Matrix<F> const &matrix);
  void invert(CTF::Matrix<F> &pseudoInverse, double accuracy = 1e-3);

  static void test(CTF::World *world);
protected:
  CTF::Matrix<F> matrix, inverse, square;
};

#endif

