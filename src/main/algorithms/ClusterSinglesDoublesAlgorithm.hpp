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

#ifndef CLUSTER_SINGLES_DOUBLES_ALGORITHM_DEFINED 
#define CLUSTER_SINGLES_DOUBLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/TensorUnion.hpp>
#include <util/SharedPointer.hpp>

#include <string>
#include <initializer_list>

namespace cc4s {
  /**
   * \brief Contains all the necessary tools for an algorithm with
   * singles and doubles amplitudes. It calculates the energy from the amplitudes
   * \f$T_{a}^{i}\f$ and \f$T_{ab}^{ij}\f$ and the Coulomb integrals \f$V_{ij}^{ab}\f$. For
   * calculating the amplitudes it calls the iteration routine of the actual algorithm.
   **/
  class ClusterSinglesDoublesAlgorithm: public Algorithm {
  public:
    /**
     * \brief Calculates the energy of a ClusterSinglesDoubles algorithm
     */
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;

    /**
     * \brief Returns the abbreviation of the concrete algorithm, e.g. "Ccd",
     * "Dcd".
     */
    virtual std::string getAbbreviation() = 0;

    /**
     * \brief Defines the default number of iterations (16).
     */
    static int constexpr DEFAULT_MAX_ITERATIONS = 16;
    static Real<64> constexpr DEFAULT_ENERGY_CONVERGENCE = 1E-7;
    static Real<64> constexpr DEFAULT_AMPLITUDES_CONVERGENCE = 1E-6;

    static Real<64> constexpr DEFAULT_LEVEL_SHIFT = 0.0;

  protected:
    Ptr<MapNode> arguments, energy;

    template <typename F, typename TE>
    Ptr<MapNode> run();

    /**
     * \brief Computes and returns the residuum of the given amplitudes.
     **/
    virtual Ptr<TensorUnion<Real<>, DefaultDryTensorEngine>>getResiduum(
      const int iteration,
      const Ptr<TensorUnion<Real<>,DefaultDryTensorEngine>> &amplitudes
    ) = 0;

    /**
     * \brief Computes and returns the residuum of the given amplitudes.
     **/
    virtual Ptr<TensorUnion<Complex<>, DefaultDryTensorEngine>>getResiduum(
      const int iteration,
      const Ptr<TensorUnion<Complex<>, DefaultDryTensorEngine>> &amplitudes
    ) = 0;
    /**
     * \brief Computes and returns the residuum of the given amplitudes.
     **/
    virtual Ptr<TensorUnion<Real<>, DefaultTensorEngine>>
    getResiduum(
      const int iteration,
      const Ptr<TensorUnion<Real<>, DefaultTensorEngine>> &amplitudes
    ) = 0;

    /**
     * \brief Computes and returns the residuum of the given amplitudes.
     **/
    virtual Ptr<TensorUnion<Complex<>, DefaultTensorEngine>>
    getResiduum(
      const int iteration,
      const Ptr<TensorUnion<Complex<>, DefaultTensorEngine>> &amplitudes
    ) = 0;

    /**
     * \brief Computes and returns the energy of the given amplitudes.
     **/
    template <typename F, typename TE>
    F getEnergy( const Ptr<TensorUnion<F,TE>> &amplitdues
               , const bool finalReport = false
               );

    /**
     * \brief Calculates an improved estimate of the amplitudes provided
     * the given the residuum.
     * \f$T_{ij\ldots}^{ab\ldots} = \frac{R_{ij\ldots}^{a\ldots}}
       {-\Delta_{ij\ldots}^{ab\ldots}}\f$
     * with \f$\Delta_{ij\ldots}^{ab\ldots} =
       \varepsilon_i+\ldots-\varepsilon_a-\ldots\f$.
     * \param[inout] residuum Fock vector, overwritten with new amplitudes.
     * \param[in] amplitudes Fock vector, previous amplitudes
     **/
    template <typename F, typename TE>
    void estimateAmplitudesFromResiduum(
      const Ptr<TensorUnion<F,TE>> &residuum,
      const Ptr<TensorUnion<F,TE>> &amplitudes
    );

    /**
     * \brief Calculates eps_a+eps_b+...-eps_i-eps_j-... into D^ab..._ij...
     **/
    template <typename F, typename TE>
    Ptr<Tensor<F,TE>> calculateExcitationEnergies(
      const std::vector<size_t> &lens, const std::string &indices
    );

    template <typename F, typename TE>
    Ptr<TensorUnion<F,TE>> createAmplitudes(
      std::initializer_list<std::initializer_list<size_t>> amplitudeLens,
      std::initializer_list<std::string> amplitudeIndices
    );

    /**
     * \brief The abbreviation of the algorithm in capital letters.
     **/
    std::string getCapitalizedAbbreviation();

    std::string getDataName(const std::string &type, const std::string &data);

    bool restart = false;
  };
}

#endif

