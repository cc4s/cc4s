#ifndef RANDOM_TENSOR_DEFINED
#define RANDOM_TENSOR_DEFINED

#include <math/Complex.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <complex>
#include <random>

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

  class DefaultRandomEngine: public std::mt19937 {
  public:
    DefaultRandomEngine() {
      seed(Cc4s::world->rank);
    }
  };

/*
    std::normal_distribution<double> normalDistribution(0.0, 1.0);
*/
  template <
    typename F, typename Distribution, typename RandomEngine
  >
  void setRandomTensor(
    CTF::Tensor<F> &t, Distribution &distribution, RandomEngine &randomEngine
  ) {
    int64_t indicesCount, *indices;
    F *values;
    t.read_local(&indicesCount, &indices, &values);
    for (int64_t i(0); i < indicesCount; ++i) {
      setRandom(values[i], distribution, randomEngine);
    }
    t.write(indicesCount, indices, values);
    free(indices); free(values);
  }
}

#endif

