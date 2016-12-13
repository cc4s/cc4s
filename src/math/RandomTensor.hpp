#ifndef RANDOM_TENSOR_DEFINED
#define RANDOM_TENSOR_DEFINED

#include <math/Complex.hpp>
#include <ctf.hpp>
#include <complex>

namespace cc4s {
  template <
    typename Distribution, typename RandomEngine
  >
  inline void setRandom(
    double &value, Distribution &distribution, RandomEngine &randomEngine
  ) {
//#ifdef INTEL_COMPILER
//    value = distribution(randomEngine);
//    value = -1.0 + 2.0*rand() / RAND_MAX; // distribution(randomEngine);
//#else
    value = distribution(randomEngine);
//#endif
  }

  template <
    typename Distribution, typename RandomEngine
  >
  inline void setRandom(
    complex &value, Distribution &distribution, RandomEngine &randomEngine
  ) {
//#ifdef INTEL_COMPILER
//    value.real() = -1.0 + 2.0*rand() / RAND_MAX; // distribution(randomEngine);
//    value.imag() = -1.0 + 2.0*rand() / RAND_MAX; // distribution(randomEngine);
//#else
    value.real(distribution(randomEngine));
    value.imag(distribution(randomEngine));
//#endif
  }

  template <typename F>
  void setRandomTensor(CTF::Tensor<F> &t);
}

#endif

