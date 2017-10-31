/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CLUSTER_SINGLES_DOUBLES_ALGORITHM_DEFINED 
#define CLUSTER_SINGLES_DOUBLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/FockVector.hpp>
#include <tcc/DryTensor.hpp>
#include <util/SharedPointer.hpp>

#include <ctf.hpp>

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
    ClusterSinglesDoublesAlgorithm(
      std::vector<Argument> const &argumentList
    );
    virtual ~ClusterSinglesDoublesAlgorithm();
    /**
     * \brief Calculates the energy of a ClusterSinglesDoubles algorithm
     */
    virtual void run();

    // TODO: dryRun

    /**
     * \brief Returns the abbreviation of the concrete algorithm, e.g. "Ccd",
     * "Dcd".
     */
    virtual std::string getAbbreviation() = 0;

    /**
     * \brief Defines the default number of iterations (16).
     */
    static int constexpr DEFAULT_MAX_ITERATIONS = 16;

    static double constexpr DEFAULT_LEVEL_SHIFT = 0.0;

  protected:
    template <typename F>
    F run();

    /**
     * \brief Computes and returns the residuum of the given amplitudes
     **/
    virtual PTR(FockVector<double>) getResiduum(
      const int iteration, const PTR(const FockVector<double>) &amplitudes
    ) = 0;

    /**
     * \brief Computes and returns the residuum of the given amplitudes
     **/
    virtual PTR(FockVector<complex>) getResiduum(
      const int iteration, const PTR(const FockVector<complex>) &amplitudes
    ) = 0;

    /**
     * \brief Computes and returns the energy of the given amplitudes.
     **/
    template <typename F>
    F getEnergy(const PTR(const FockVector<F>) &amplitdues);

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
    template <typename F>
    void estimateAmplitudesFromResiduum(
      const PTR(FockVector<F>) &residuum,
      const PTR(const FockVector<F>) &amplitudes
    );

    /**
     * \brief Calculates eps_a+eps_b+...-eps_i-eps_j-... into D^ab..._ij...
     **/
    template <typename F>
    void calculateExcitationEnergies(
      CTF::Tensor<F> &D, const std::string &indices
    );

    /**
     * \brief Dry run for amplitudesFromResiduum.
     * \param[in] R residuum tensor.
     **/
    template <typename F>
    void dryAmplitudesFromResiduum(cc4s::DryTensor<F> &R);

    template <typename F>
    PTR(FockVector<F>) createAmplitudes(
      std::initializer_list<std::string> amplitudeNames,
      std::initializer_list<std::initializer_list<int>> amplitudeLens,
      std::initializer_list<std::string> amplitudeIndices
    );

    template <typename F>
    void storeAmplitudes(
      const PTR(const FockVector<F>) &amplitudes,
      std::initializer_list<std::string> names
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
    CTF::Tensor<double> *sliceCoupledCoulombIntegrals(
      const PTR(const FockVector<double>) &amplitudes,
      int a, int b, int integralsSliceSize
    );

    CTF::Tensor<complex> *sliceCoupledCoulombIntegrals(
      const PTR(const FockVector<complex>) &amplitudes,
      int a, int b, int integralsSliceSize
    );


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
    CTF::Tensor<double> *sliceAmplitudesFromCoupledCoulombFactors(
      const PTR(const FockVector<double>) &amplitudes,
      int a, int b, int factorsSliceSize
    );
    CTF::Tensor<complex> *sliceAmplitudesFromCoupledCoulombFactors(
      const PTR(const FockVector<complex>) &amplitudes,
      int a, int b, int factorsSliceSize
    );

    /**
     * \brief Adds the given slice of the residuum tensor Rxyij to the
     * entire residuum tensor Rabij at the respective index range.
     * \param[in] a0 1st sliced dimension (x).
     * \param[in] b0 2nd sliced dimension (y).
     * \param[in] Rxyij sliced residuum
     * \param[in] Rabij entire residuum.
     */
    template <typename F>
    void sliceIntoResiduum(
      CTF::Tensor<F> &Rxyij, int a0, int b0, CTF::Tensor<F> &Rabij
    );


    /**
     * \brief The abbreviation of the algorithm in capital letters.
     **/
    std::string getCapitalizedAbbreviation();

    std::string getDataName(const std::string &type, const std::string &data);
  };
}

#endif

