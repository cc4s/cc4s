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

#ifndef COUPLED_CLUSTER_SINGLES_DOUBLES_DEFINED
#define COUPLED_CLUSTER_SINGLES_DOUBLES_DEFINED


#include <algorithms/CoupledClusterSinglesDoublesReference.hpp>

#include <util/SharedPointer.hpp>

namespace cc4s {
  // this algorithm is now based on the ClusterSinglesDoublesAlgorithm
  // inheriting its iteration and slicing functionality.
  // Only the abstract (left out) methods getAbbreviation and iterate have
  // to be implemented.
  /**
   * \brief Implements the iteration routine for the Drccd method. Calculates the
   * amplitudes \f$T_{ab}^{ij}\f$ from the Coulomb Integrals \f$V_{ij}^{ab}\f$
   * in a \f$ \mathcal{O}(N^{6}) \f$ implementation.
   */
  class CoupledClusterSinglesDoubles:
    public CoupledClusterSinglesDoublesReference {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CoupledClusterSinglesDoubles)
    /**
     * \brief Returns the abbreviation of the routine (CCSD).
     * \return abbreviation of the routine
     */
    std::string getAbbreviation() override { return "Ccsd"; }
  protected:
    /**
     * \brief Implements the iterate method with the CCSD iteration.
     * \param[in] i Iteration number
     */
    Ptr<TensorUnion<Real<>, DefaultDryTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<TensorUnion<Real<>, DefaultDryTensorEngine>> &amplitudes
    ) override;
    Ptr<TensorUnion<Complex<>, DefaultDryTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<TensorUnion<Complex<>, DefaultDryTensorEngine>> &amplitudes
    ) override;
    Ptr<TensorUnion<Real<>, DefaultTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<TensorUnion<Real<>, DefaultTensorEngine>> &amplitudes
    ) override;
    Ptr<TensorUnion<Complex<>, DefaultTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<TensorUnion<Complex<>, DefaultTensorEngine>> &amplitudes
    ) override;

    template <typename TE>
    Ptr<TensorUnion<Real<>,TE>> getResiduum(
      const int iteration, const Ptr<TensorUnion<Real<>,TE>> &amplitudes
    );
    template <typename TE>
    Ptr<TensorUnion<Complex<>,TE>> getResiduum(
      const int iteration, const Ptr<TensorUnion<Complex<>,TE>> &amplitudes
    );
  };
}

#endif

