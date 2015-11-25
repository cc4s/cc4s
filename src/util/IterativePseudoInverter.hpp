#ifndef ITERATIVE_PSEUDO_INVERTER_DEFINED
#define ITERATIVE_PSEUDO_INVERTER_DEFINED

#include <util/Complex.hpp>
#include <ctf.hpp>
#include <random>

template <typename F>
class IterativePseudoInverter {
public:
  IterativePseudoInverter(CTF::Matrix<F> const &matrix);
  CTF::Matrix<F> &invert(double accuracy = 1e-10);

  static void test(CTF::World *world);
protected:
  void iterate(double accuracy = 1e-10);
  void iterateQuadratically(double accuracy = 1e-10);
  static void setRandom(
    F &value,
    std::mt19937 &random, std::normal_distribution<double> &normalDistribution
  );
  static void generateHilbertMatrix(CTF::Matrix<F> &matrix);

  CTF::Matrix<F> matrix, square, inverse;
  double alpha;
};

#endif
