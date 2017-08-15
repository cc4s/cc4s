/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CLUSTER_SINGLES_DOUBLES_ALGORITHM_DEFINED 
#define CLUSTER_SINGLES_DOUBLES_ALGORITHM_DEFINED

#include <algorithms/ClusterDoublesAlgorithm.hpp>
#include <mixers/Mixer.hpp>
#include <string>
#include <ctf.hpp>
#include <tcc/DryTensor.hpp>

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

  protected:
    template <typename F>
    F run();

    /**
     * \brief Performs one iteration of the concrete algorithm.
     **/
    virtual void iterate(
      int i, Mixer<double> *TaiMixer, Mixer<double> *TabijMixer
    ) = 0;

    /**
     * \brief Performs one iteration of the concrete algorithm.
     **/
    virtual void iterate(
      int i, Mixer<complex> *TaiMixer, Mixer<complex> *TabijMixer
    ) = 0;

    /**
     * \brief Calculates the energy from the amplitudes currently contained
     * in the mixers. Overrides ClusterDoublesAlgorithm method.
     **/
    template <typename F>
    F calculateEnergy(Mixer<F> *TaiMixer, Mixer<F> *TabijMixer);

    /**
     * \brief Calculates the amplitudes from the current residuum and
     * returns them in-place.
     * Usually this is done by calculating
     * \f$T_{ij\ldots}^{ab\ldots} = R_{ij\ldots}^{a\ldots} / (\varepsilon_i+\ldots-\varepsilon_a-\ldots)\f$,
     * but other methods, such as level shifting may be used.
     * \param[in] R residuum tensor.
     * \param[in] indices indices into residuum tensor, e.g. "abij".
     **/
    template <typename F>
    void amplitudesFromResiduum(
      CTF::Tensor<F> &R, const std::string &indices
    );

    /**
     * \brief Dry run for amplitudesFromResiduum.
     * \param[in] R residuum tensor.
     **/
    template <typename F>
    void dryAmplitudesFromResiduum(cc4s::DryTensor<F> &R);

    template <typename F>
    Mixer<F> *createMixer(const std::string &type, std::vector<int> shape);

    template <typename F>
    void storeAmplitudes(Mixer<F> *mixer, const std::string &type);

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
      Mixer<double> *TaiMixer, int a, int b, int integralsSliceSize
    );

    CTF::Tensor<complex> *sliceCoupledCoulombIntegrals(
      Mixer<complex> *TaiMixer, int a, int b, int integralsSliceSize
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
      Mixer<double> *TaiMixer, Mixer<double> *TabijMixer,
      int a, int b, int factorsSliceSize
    );
    CTF::Tensor<complex> *sliceAmplitudesFromCoupledCoulombFactors(
      Mixer<complex> *TaiMixer, Mixer<complex> *TabijMixer,
      int a, int b, int factorsSliceSize
    );

    /**
     * \brief Calculates and returns tensor Fabij of the residuum
     * from the dressed Coulomb factors. 
     * The caller is responsible for deleting the dynamically allocated
     * result tensor. 
     * \param[out] Fabij sliced Residuum
     */
    CTF::Tensor<double> *amplitudesFromCoupledCoulombFactors(
      Mixer<double> *TaiMixer, Mixer<double> *TabijMixer
    );
    CTF::Tensor<complex> *amplitudesFromCoupledCoulombFactors(
      Mixer<complex> *TaiMixer, Mixer<complex> *TabijMixer
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

