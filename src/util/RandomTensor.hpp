#ifndef RANDOM_TENSOR_DEFINED
#define RANDOM_TENSOR_DEFINED

#include <util/Complex.hpp>
#include <ctf.hpp>

namespace cc4s {
  template <
    typename Distribution, typename RandomEngine
  >
  void setRandom(
    double &value, Distribution &distribution, RandomEngine &randomEngine
  ) {
    value = distribution(randomEngine);
  }

  template <
    typename Distribution, typename RandomEngine
  >
  void setRandom(
    complex &value, Distribution &distribution, RandomEngine &randomEngine
  ) {
    value.real(distribution(randomEngine));
    value.imag(distribution(randomEngine));
  }

  template <typename F>
  void setRandomTensor(CTF::Tensor<F> &t);
}

#endif

