/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef DISTRIBUTED_SAMPLED_VARIABLE_DEFINED
#define DISTRIBUTED_SAMPLED_VARIABLE_DEFINED

#include <math/SampledVariable.hpp>
#include <util/MpiCommunicator.hpp>

namespace cc4s {
  template <typename F=Real<>>
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

