/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CLUSTER_SINGLES_DOUBLES_ALGORITHM_DEFINED 
#define CLUSTER_SINGLES_DOUBLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
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
    virtual Ptr<FockVector<Real<>, DryTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<const FockVector<Real<>, DryTensorEngine>> &amplitudes
    ) = 0;

    /**
     * \brief Computes and returns the residuum of the given amplitudes.
     **/
    virtual Ptr<FockVector<Complex<>, DryTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<const FockVector<Complex<>, DryTensorEngine>> &amplitudes
    ) = 0;
    /**
     * \brief Computes and returns the residuum of the given amplitudes.
     **/
    virtual Ptr<FockVector<Real<>, DefaultTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<const FockVector<Real<>, DefaultTensorEngine>> &amplitudes
    ) = 0;

    /**
     * \brief Computes and returns the residuum of the given amplitudes.
     **/
    virtual Ptr<FockVector<Complex<>, DefaultTensorEngine>> getResiduum(
      const int iteration,
      const Ptr<const FockVector<Complex<>, DefaultTensorEngine>> &amplitudes
    ) = 0;

    /**
     * \brief Computes and returns the energy of the given amplitudes.
     **/
    template <typename F, typename TE>
    F getEnergy( const Ptr<const FockVector<F,TE>> &amplitdues
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
      const Ptr<FockVector<F,TE>> &residuum,
      const Ptr<const FockVector<F,TE>> &amplitudes
    );

    /**
     * \brief Calculates eps_a+eps_b+...-eps_i-eps_j-... into D^ab..._ij...
     **/
    template <typename F, typename TE>
    Ptr<Tensor<F,TE>> calculateExcitationEnergies(
      const std::vector<size_t> &lens, const std::string &indices
    );

    template <typename F, typename TE>
    Ptr<FockVector<F,TE>> createAmplitudes(
      std::initializer_list<std::string> amplitudeComponent,
      std::initializer_list<std::initializer_list<size_t>> amplitudeLens,
      std::initializer_list<std::string> amplitudeIndices
    );

    template <typename F, typename TE>
    Ptr<MapNode> storeAmplitudes(
      const Ptr<MapNode> &arguments,
      const Ptr<const FockVector<F,TE>> &amplitudes
    );
    template <typename F, typename TE>
    Ptr<MapNode> storeAmplitudesComponent(
      const Ptr<Tensor<F,TE>> &component
    );

    /**
     * \brief Calculates and returns one slice Xxycd of the Coulomb integrals \f$V_{cd}^{ab}\f$
     * coupled to the singles amplitudes.
     * The indices x and y are restricted to the
     * range {No+a, ..., No+a+No-1} and {No+b, ..., No+b+No-1}, respectively.
     * The caller is responsible for deleting the dynamically allocated
     * result tensor. 
     * \param[in] a 1st sliced dimension (x).
     * \param[in] b 2nd sliced dimension (y).
     * \param[in] integralsSliceSize slicing rank.
     * \param[out] Xxycd sliced coupled Coulomb integrals Xabcd
     */
/*
    CTF::Tensor<double> *sliceCoupledCoulombIntegrals(
      const Ptr<const FockVector<double>> &amplitudes,
      int a, int b, int integralsSliceSize
    );

    CTF::Tensor<complex> *sliceCoupledCoulombIntegrals(
      const Ptr<const FockVector<complex>> &amplitudes,
      int a, int b, int integralsSliceSize
    );
*/


    /**
     * \brief Calculates and returns one slice Fabij of the residuum
     * from the dressed Coulomb factors. The slice is computed from
     * Rx and Ry and are restricted to the
     * range {a, ..., factorsSliceSize+a-1} and {b, ..., factorsSliceSize+b-1}, respectively.
     * The caller is responsible for deleting the dynamically allocated
     * result tensor. 
     * \param[in] a 1st sliced dimension (Rx).
     * \param[in] b 2nd sliced dimension (Ry).
     * \param[in] factorsSliceSize slicing rank of NR.
     * \param[out] Fabij sliced Residuum
     */
/*
    CTF::Tensor<double> *sliceAmplitudesFromCoupledCoulombFactors(
      const Ptr<const FockVector<double>> &amplitudes,
      int a, int b, int factorsSliceSize
    );
    CTF::Tensor<complex> *sliceAmplitudesFromCoupledCoulombFactors(
      const Ptr<const FockVector<complex>> &amplitudes,
      int a, int b, int factorsSliceSize
    );
*/

    /**
     * \brief Adds the given slice of the residuum tensor Rxyij to the
     * entire residuum tensor Rabij at the respective index range.
     * \param[in] a0 1st sliced dimension (x).
     * \param[in] b0 2nd sliced dimension (y).
     * \param[in] Rxyij sliced residuum
     * \param[in] Rabij entire residuum.
     */
/*
    template <typename F>
    void sliceIntoResiduum(
      CTF::Tensor<F> &Rxyij, int a0, int b0, CTF::Tensor<F> &Rabij
    );
*/

    /**
     * \brief The abbreviation of the algorithm in capital letters.
     **/
    std::string getCapitalizedAbbreviation();

    std::string getDataName(const std::string &type, const std::string &data);

    bool restart = false;
  };
}

#endif

