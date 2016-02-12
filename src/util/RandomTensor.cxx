#include <util/RandomTensor.hpp>

#include <util/Complex.hpp>
#include <util/Log.hpp>
#include <random>

/*
template <
  typename Distribution,
  typename RandomEngine
>
void cc4s::setRandom(
  double &value, Distribution &distribution, RandomEngine &randomEngine
) {
  value = distribution(randomEngine);
}

template <
  typename Distribution,
  typename RandomEngine
>
void cc4s::setRandom(
  complex &value, Distribution &distribution, RandomEngine &randomEngine
) {
  value.real(distribution(randomEngine));
  value.imag(distribution(randomEngine));
}
*/

template <typename F>
void cc4s::setRandomTensor(CTF::Tensor<F> &t) {
  int64_t indicesCount, *indices;
  F *values;
  std::mt19937 random;
  random.seed(t.wrld->rank);
  std::normal_distribution<double> normalDistribution(0.0, 1.0);
  t.read_local(&indicesCount, &indices, &values);
  for (int64_t i(0); i < indicesCount; ++i) {
    setRandom(values[i], normalDistribution, random);
  }
  t.write(indicesCount, indices, values);
  free(indices); free(values);
}

// instantiate
template
void cc4s::setRandomTensor(CTF::Tensor<double> &t);
template
void cc4s::setRandomTensor(CTF::Tensor<complex> &t);

