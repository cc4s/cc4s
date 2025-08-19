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

#ifndef CCSD_DEFINED
#define CCSD_DEFINED

#include <algorithms/coupledcluster/CoupledClusterMethod.hpp>
#include <SharedPointer.hpp>

namespace cc4s {
  /**
   * \brief Implements the iteration routine for the ccsd method. Calculates the
   * amplitudes \f$T_{ab}^{ij}\f$ from the Coulomb Integrals \f$V_{ij}^{ab}\f$
   * in a \f$ \mathcal{O}(N^{6}) \f$ implementation.
   */

  template <typename F, typename TE>
  class Ccsd : public CoupledClusterMethod<F, TE> {
  public:
    Ccsd(const Ptr<MapNode> &arguments)
      : CoupledClusterMethod<F, TE>(arguments) {}

    std::string getName() override { return "Ccsd"; }

    static CoupledClusterMethodRegistrar<F, TE, Ccsd<F, TE>> registrar_;

    std::string describeOptions() override;

    Ptr<TensorSet<F, TE>> getResiduum(const Ptr<TensorSet<F, TE>> &amplitudes) override;
  };

}

#endif

