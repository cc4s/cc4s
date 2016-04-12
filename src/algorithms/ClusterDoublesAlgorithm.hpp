/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CLUSTER_DOUBLES_ALGORITHM_DEFINED 
#define CLUSTER_DOUBLES_ALGORITHM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <mixers/Mixer.hpp>
#include <string>
#include <ctf.hpp>

namespace cc4s {
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
     * \brief Returns the abbreviation of the concrete algorithm, e.g. "CCD",
     * "DCD"
     */
    virtual std::string getAbbreviation() = 0;

    static int64_t constexpr DEFAULT_MAX_ITERATIONS = 16;

  protected:
    /**
     * \brief The mixer for the doubles amplitudes.
     */
    Mixer<double> *TabijMixer;

    /**
     * \brief Performs one iteration of the concrete algorithm.
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
     */
    void doublesAmplitudesFromResiduum(CTF::Tensor<> &Rabij);

    /**
     * \brief Calculates and returns one slice Vxycd of the Coulomb integrals
     * from the Coulomb vertex. The indices x and y are restricted to the
     * range {No+a, ..., No+a+No-1} and {No+b, ..., No+b+No-1}, respectively.
     * The caller is responsible for deleting the dynamically allocated
     * result tensor. 
     */
    CTF::Tensor<> *sliceCoulombIntegrals(int a, int b, int sliceRank);
    /**
     * \brief Adds the given slice of the residuum tensor Rxyij to the
     * entire residuum tensor Rabij at the respective index range.
     */
    void sliceIntoResiduum(
      CTF::Tensor<> &Rxyij, int a0, int b0, CTF::Tensor<> &Rabij
    );
  };
}

#endif

