#ifndef ITERATIVE_PSEUDO_INVERTER_DEFINED
#define ITERATIVE_PSEUDO_INVERTER_DEFINED

#include <util/Complex.hpp>
#include <ctf.hpp>

class IterativePseudoInverter {
public:
  IterativePseudoInverter(CTF::Matrix<complex> const &matrix);
  void invert(CTF::Matrix<complex> &pseudoInverse, double accuracy);

protected:
  CTF::Matrix<complex> matrix, square;
};

#endif

