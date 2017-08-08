/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CLUSTER_DOUBLES_ALGORITHM_DEFINED 
#define CLUSTER_DOUBLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <mixers/Mixer.hpp>
#include <string>
#include <ctf.hpp>
#include <tcc/DryTensor.hpp>

namespace cc4s {
  /**
   * \brief Contains all the necessary tools for an algorithm with
   * only doubles amplitudes. It calculates the energy from the amplitudes
   * \f$T_{ab}^{ij}\f$ and the Coulomb integrals \f$V_{ij}^{ab}\f$. For
   * calculating the amplitudes it calls the iteration routine of the actual algorithm.
   */
  class ClusterDoublesAlgorithm: public Algorithm {
  public:
    ClusterDoublesAlgorithm(
      std::vector<Argument> const &argumentList
    );
    virtual ~ClusterDoublesAlgorithm();
    /**
     * \brief Calculates the energy of a ClusterDoubles algorithm
     */
    virtual void run();

    /**
     * \brief Performs a Dry Run
     */
    virtual void dryRun();
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
    /**
     * \brief Evaluates and returns the energy according to the respective
     * ClusterDoubles algorithm. An instantiation of this method
     * is called from the run method with typename F being either double
     * or complex, depending on which type of integrals are present.
     **/
    template <typename F>
    F run();

    /**
     * \brief Evaluates a dry run and returns 0 according to the respective
     * ClusterDoubles algorithm. An instantiation of this method
     * is called from the dryRun method with typename F being either double
     * or complex, depending on which type of integrals are present.
     **/
    template <typename F>
    F dryRun();

    /**
     * \brief Performs one iteration of the concrete algorithm
     * when operating with real valued Coulomb integrals.
     * \param[in] i Iteration number
     * \param[in] TaiMixer Mixer for the singles amplitudes, may be null
     * \param[in] TabijMixer Mixer for the doubles amplitudes.
     */
    virtual void iterate(
      int i, Mixer<double> *TaiMixer, Mixer<double> *TabijMixer
    ) = 0;
    /**
     * \brief Performs one iteration of the concrete algorithm
     * when operating with complex valued Coulomb integrals.
     * \param[in] i Iteration number
     * \param[in] TaiMixer Mixer for the singles amplitudes, may be null
     * \param[in] TabijMixer Mixer for the doubles amplitudes.
     */
    virtual void iterate(
      int i, Mixer<complex> *TaiMixer, Mixer<complex> *TabijMixer
    ) = 0;

    /**
     * \brief Performs a dry run of one iteration of the concrete algorithm
     * when operating with real valued Coulomb integrals.
     * The base class does not perform accounting and writes a warning about
     * that.
     * Note that there is no DryMixer currently.
     * \param[in] TaiMixer singles amplitudes, may be null
     * \param[in] TabijMixer doubles amplitudes.
     */
    virtual void dryIterate(
      DryTensor<double> *TaiMixer, DryTensor<double> *TabijMixer
    );
    /**
     * \brief Performs a dry run of one iteration of the concrete algorithm
     * when operating with complex valued Coulomb integrals.
     * The base class does not perform accounting and writes a warning about
     * that.
     * \param[in] TaiMixer Mixer for the singles amplitudes, may be null
     * \param[in] TabijMixer Mixer for the doubles amplitudes.
     */
    virtual void dryIterate(
      DryTensor<complex> *TaiMixer, DryTensor<complex> *TabijMixer
    );

    template <typename F>
    Mixer<F> *createDoublesMixer();

    template <typename F>
    void storeDoublesAmplitudes(Mixer<F> *TabijMixer);

    /**
     * \brief Calculates the energy from the amplitudes currently contained
     * in the mixer.
     * \param[in] TabijMixer doubles amplitudes mixer.
     **/
    template <typename F>
    F calculateEnergy(Mixer<F> *TabijMixer);

    /**
     * \brief Calculates the doubles amplitudes from the given doubles residuum
     * and returns them in place of the residuum tensor.
     * Usually this is done by calculating
     * \f$T_{ij}^{ab} = R_{ij}^{ab} / (\varepsilon_i+\varepsilon_j-\varepsilon_a-\varepsilon_b)\f$,
     * but other methods, such as level shifting may be used.
     * \param[in] Rabij doubles amplitudes mixer.
     **/
    template <typename F>
    void doublesAmplitudesFromResiduum(CTF::Tensor<F> &Rabij);

    /**
     * \brief Dry run for doublesAmplitudesFromResiduum.
     * \param[in] Rabij Residuum Tensor.
     */
    void dryDoublesAmplitudesFromResiduum(cc4s::DryTensor<> &Rabij);

    /**
     * \brief The abbreviation of the algorithm in capital letters.
     **/
    std::string getCapitalizedAbbreviation();

    /**
     * \brief Constructs the doubles amplitudes name from the abbreviation.
     **/
    std::string getDoublesAmplitudesName();
  };
}

#endif

