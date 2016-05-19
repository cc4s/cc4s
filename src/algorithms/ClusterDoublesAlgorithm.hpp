/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CLUSTER_DOUBLES_ALGORITHM_DEFINED 
#define CLUSTER_DOUBLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <mixers/Mixer.hpp>
#include <string>
#include <ctf.hpp>

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

    /** \brief The occupied orbital energies  */
    CTF::Tensor<> *epsi;
    /** \brief The virtual orbital energies  */
    CTF::Tensor<> *epsa;
    /** \brief The Coulomb integrals Vabij  */
    CTF::Tensor<> *Vabij;
    /** \brief The Coulomb Vertex GammaGpq  */
    CTF::Tensor<complex> *GammaGpq;

    /**
     * \brief Performs a Dry Run
     */
    virtual void dryRun();
    /**
     * \brief Returns the abbreviation of the concrete algorithm, e.g. "CCD",
     * "DCD"
     */
    virtual std::string getAbbreviation() = 0;

    /**
     * \brief Defines the default number of iterations (16).
     */
    static int constexpr DEFAULT_MAX_ITERATIONS = 16;

  protected:
    /**
     * \brief The mixer for the doubles amplitudes.
     */
    Mixer<double> *TabijMixer;

    /**
     * \brief Performs one iteration of the concrete algorithm.
     * \param[in] i Iteration number
     */
    virtual void iterate(int i) = 0;

    /**
     * \brief Performs a dry run of one iteration of the concrete algorithm.
     * The base class does not perform accounting and writes a warning about
     * that.
     */
    virtual void dryIterate();

    /**
     * \brief Calculates the amplitudes from the current residuum and
     * returns them in-place.
     * Usually this is done by calculating
     * \f$T_{ij}^{ab} = R_{ij}^{ab} / (\varepsilon_i+\varepsilon_j-\varepsilon_b-\varepsilon_b)\f$,
     * but other methods, such as level shifting may be used.
     * \param[in] Rabij Residuum Tensor.
     */
    void doublesAmplitudesFromResiduum(CTF::Tensor<> &Rabij);

    /**
     * \brief Calculates and returns one slice Vxycd of the Coulomb integrals
     * from the Coulomb vertex. The indices x and y are restricted to the
     * range {No+a, ..., No+a+No-1} and {No+b, ..., No+b+No-1}, respectively.
     * The caller is responsible for deleting the dynamically allocated
     * result tensor. 
     * \param[in] a 1st sliced dimension (x).
     * \param[in] b 2nd sliced dimension (y).
     * \param[in] sliceRank slicing rank.
     * \param[out] Vxycd sliced Coulomb integrals Vabcd
     */
    CTF::Tensor<> *sliceCoulombIntegrals(int a, int b, int sliceRank);

    /**
     * \brief Adds the given slice of the residuum tensor Rxyij to the
     * entire residuum tensor Rabij at the respective index range.
     * \param[in] a0 1st sliced dimension (x).
     * \param[in] b0 2nd sliced dimension (y).
     * \param[in] Rxyij sliced residuum
     * \param[in] Rabij entire residuum.
     */
    void sliceIntoResiduum(
      CTF::Tensor<> &Rxyij, int a0, int b0, CTF::Tensor<> &Rabij
    );
  };
}

#endif

