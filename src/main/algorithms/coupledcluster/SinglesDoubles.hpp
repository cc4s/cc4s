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

#ifndef SINGLES_DOUBLES_DEFINED
#define SINGLES_DOUBLES_DEFINED


#include <algorithms/coupledcluster/CoupledClusterMethod.hpp>

#include <util/SharedPointer.hpp>

namespace cc4s {
  /**
   * \brief Implements the iteration routine for the ccsd method. Calculates the
   * amplitudes \f$T_{ab}^{ij}\f$ from the Coulomb Integrals \f$V_{ij}^{ab}\f$
   * in a \f$ \mathcal{O}(N^{6}) \f$ implementation.
   */
  template <typename F, typename TE>
  class SinglesDoubles {};

  template <typename TE>
  class SinglesDoubles<Real<>,TE>: public CoupledClusterMethod<Real<>,TE> {
  public:
    SinglesDoubles(
      const Ptr<MapNode> &arguments
    ): CoupledClusterMethod<Real<>,TE>(arguments) {
    }
    std::string getName() override { return "SinglesDoubles"; } \
    static CoupledClusterMethodRegistrar<
      Real<>,TE,SinglesDoubles<Real<>,TE>
    > registrar_;

    /**
     * \brief Implements the iterate method with the CCSD iteration.
     * \param[in] i Iteration number
     */
    Ptr<TensorUnion<Real<>,TE>> getResiduum(
      const int iteration, const bool restart,
      const Ptr<TensorUnion<Real<>,TE>> &amplitudes
    ) override;
  };

  template <typename TE>
  class SinglesDoubles<Complex<>,TE>: public CoupledClusterMethod<Complex<>,TE> {
  public:
    SinglesDoubles(
      const Ptr<MapNode> &arguments
    ): CoupledClusterMethod<Complex<>,TE>(arguments) {
    }
    std::string getName() override { return "SinglesDoubles"; } \
    static CoupledClusterMethodRegistrar<
      Complex<>,TE,SinglesDoubles<Complex<>,TE>
    > registrar_;

    /**
     * \brief Implements the iterate method with the CCSD iteration.
     * \param[in] i Iteration number
     */
    Ptr<TensorUnion<Complex<>,TE>> getResiduum(
      const int iteration, const bool restart,
      const Ptr<TensorUnion<Complex<>,TE>> &amplitudes
    ) override;
  };

/*
  template <typename TE>
  Ptr<TensorUnion<Real<>,TE>> SinglesDoubles<Real<>,TE>::getResiduum(
    const int iteration, const bool restart,
    const Ptr<TensorUnion<Real<>,TE>> &amplitudes
  );
*/
}

#endif

