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
     * \brief The abbreviation of the algorithm in capital letters.
     **/
    std::string getCapitalizedAbbreviation();

    std::string getDataName(const std::string &type, const std::string &data);
  };
}

#endif

