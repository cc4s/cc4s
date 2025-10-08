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

#ifndef LAMBDA_CCSD_DEFINED
#define LAMBDA_CCSD_DEFINED

#include <algorithms/coupledcluster/CoupledClusterMethod.hpp>
#include <SharedPointer.hpp>

namespace cc4s {

  template <typename F, typename TE>
  class LambdaCcsd: public CoupledClusterMethod<F,TE> {
  protected:
    // Persistent cache: keeps intermediates between iterations
    Ptr<TensorSet<F,TE>> intermediates;

  public:
    LambdaCcsd(const Ptr<MapNode> &arguments)
      : CoupledClusterMethod<F,TE>(arguments), intermediates(nullptr)
    {
    }  

  /*
  public:
    LambdaCcsd(
      const Ptr<MapNode> &arguments
    ): CoupledClusterMethod<F,TE>(arguments) {
    }
   */
    std::string getName() override { return "LambdaCcsd"; }
    static CoupledClusterMethodRegistrar<
      F,TE,LambdaCcsd<F,TE>
    > registrar_;


    Ptr<TensorSet<F,TE>> getResiduum(
      const Ptr<TensorSet<F,TE>> &amplitudes
    ) override;

    Ptr<TensorSet<F,TE>> getIntermediates(
      const Ptr<TensorSet<F,TE>> &amplitudes
    ) ;
  };
}

#endif

