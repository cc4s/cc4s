#ifndef DISTRIBUTED_SAMPLED_VARIABLE_DEFINED
#define DISTRIBUTED_SAMPLED_VARIABLE_DEFINED

#include <math/SampledVariable.hpp>
#include <util/MpiCommunicator.hpp>

namespace cc4s {
  template <typename F=double>
  class DistributedSampledVariable: public SampledVariable<F> {
  public:
    DistributedSampledVariable(
      SampledVariable<F> *globalSampledVariable_,
      MpiCommunicator *communicator_
    ):
      SampledVariable<F>(),
      globalSampledVariable(globalSampledVariable_),
      communicator(communicator_)
    {
    }
    ~DistributedSampledVariable() {
      communicator->allReduce(this->s, globalSampledVariable->s);
      communicator->allReduce(this->s2, globalSampledVariable->s2);
      communicator->allReduce(this->n, globalSampledVariable->n);
    }
  protected:
    SampledVariable<F> *globalSampledVariable;
    MpiCommunicator *communicator;
  };
}

#endif

