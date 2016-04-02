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
    virtual void run();
    /**
     * \brief Returns the abbreviation of the concrete algorithm, e.g. "CCD".
     */
    virtual std::string getAbbreviation() = 0;

    static int64_t constexpr DEFAULT_MAX_ITERATIONS = 16;

  protected:
    /**
     * \brief The mixer for the doubles amplitudes.
     */
    Mixer<double> *TabijMixer;

    /**
     * \brief erforms one iteration of the concrete algorithm.
     */
    virtual void iterate(int i) = 0;

    /**
     * \brief Calculates and returns one slice Vxycd of the Coulomb integrals
     * from the Coulomb vertex. The indices x and y are restricted to the
     * range {No+a, ..., No+a+No-1} and {No+b, ..., No+b+No-1}, respectively.
     * The caller is responsible for deleting the dynamically allocated
     * result tensor. 
     */
    CTF::Tensor<> *sliceCoulombIntegrals(int a, int b);
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

