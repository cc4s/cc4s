#include <algorithms/experimental/StochasticContraction.hpp>
#include <math/RandomGenerator.hpp>
#include <math/SampledVariable.hpp>
#include <util/DistributedSampledVariable.hpp>
#include <util/MpiCommunicator.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(StochasticContraction);

StochasticContraction::StochasticContraction(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

StochasticContraction::~StochasticContraction() {
}

cc4s::complex StochasticContraction::drawUniformWeight() {
  return complex(
    rand->nextUniform() - 0.5, rand->nextUniform()-0.5
  ) * 2.0 * std::sqrt(3.0/2.0);
}

void StochasticContraction::run() {
  // initialize random generator depending on the rank of the process
  MpiCommunicator comm(*Cc4s::world);
  rand = new RandomGenerator(Cc4s::world->rank);
  SampledVariable<complex> product;

  int64_t edges(getIntegerArgument("edges"));
  int64_t matches(getIntegerArgument("matches"));
  int64_t samplesCount(getIntegerArgument("samples"));
  LOG(1, "SC") << "edges=" << edges << ", matches=" << matches << std::endl;
  {
    DistributedSampledVariable<complex> distributedProduct(&product, &comm);
    for (int64_t n(0); n < samplesCount; ++n) {
      int i;
      complex sample(1);
      for (i = 0; i < matches; ++i) {
        complex a(drawUniformWeight());
        sample *= a * std::conj(a);
      }
      for (; i < edges; ++i) { 
        complex a(drawUniformWeight());
        complex b(drawUniformWeight());
        sample *= a * b;
      }
      distributedProduct.addSample(sample);
    }
  }

  LOG(1, "SC") << "E(product)=" << product.getMean() << std::endl;
  LOG(1, "SC") << "95% confidence interval=" << product.getMeanStdDeviation() << std::endl;
  LOG(1, "SC") << "Var(product)=" << product.getVariance() << std::endl;
  delete rand;
}

