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

#ifndef CCSDT_REFERENCE_DEFINED
#define CCSDT_REFERENCE_DEFINED

#include <algorithms/coupledcluster/CoupledClusterMethod.hpp>
#include <SharedPointer.hpp>

namespace cc4s {
  // this algorithm is now based on the ClusterSinglesDoublesAlgorithm
  // inheriting its iteration and slicing functionality.
  // Only the abstract (left out) methods getAbbreviation and iterate have
  // to be implemented.
  /**
   * \brief Implements the iteration routine for the ccsdt method.
   */
  template <typename F, typename TE>
  class Ccsdt: public CoupledClusterMethod<F,TE> {
  public:
    Ccsdt(
      const Ptr<MapNode> &arguments
    ): CoupledClusterMethod<F,TE>(arguments) {
    }
    std::string getName() override { return "Ccsdt"; }
    static CoupledClusterMethodRegistrar<
      F,TE,Ccsdt<F,TE>
    > registrar_;

    /**
     * \brief Implements the iterate method with the CCSDT iteration.
     * \param[in] i Iteration number
     */
    Ptr<TensorSet<F,TE>> getResiduum(
      const Ptr<TensorSet<F,TE>> &amplitudes
    ) override;
  };
}

#endif

