#ifndef ITERATIVE_PSEUDO_INVERSE_DEFINED
#define ITERATIVE_PSEUDO_INVERSE_DEFINED

#include <math/Complex.hpp>
#include <tcc/DryTensor.hpp>
#include <ctf.hpp>
#include <random>

namespace cc4s {
  template <typename F>
  class IterativePseudoInverse {
  public:
    IterativePseudoInverse(CTF::Tensor<F> const &matrix, double accuracy=1e-10);
    CTF::Matrix<F> &get();

    static void test(CTF::World *world);
  protected:
    void iterate(double accuracy = 1e-10);
    void iterateQuadratically(double accuracy = 1e-10);
    static void setRandom(
      F &value,
      std::mt19937 &random, std::normal_distribution<double> &normalDistribution
    );
    static void generateHilbertMatrix(CTF::Tensor<F> &matrix);

    CTF::Tensor<F> matrix;
    CTF::Matrix<F> square, inverse;
    double alpha;
  };

  template <typename F>
  class DryIterativePseudoInverse {
  public:
    DryIterativePseudoInverse(DryTensor<F> const &matrix);
    DryMatrix<F> &get();

  protected:
    DryTensor<F> matrix;
    DryMatrix<F> square, inverse;
  };
}

#endif

