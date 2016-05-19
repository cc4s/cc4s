/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CLUSTER_SINGLES_DOUBLES_ALGORITHM_DEFINED 
#define CLUSTER_SINGLES_DOUBLES_ALGORITHM_DEFINED

#include <algorithms/ClusterDoublesAlgorithm.hpp>
#include <mixers/Mixer.hpp>
#include <string>
#include <ctf.hpp>

namespace cc4s {
  /**
   * \brief Contains all the necessary tools for an algorithm with
   * singles and doubles amplitudes. It calculates the energy from the amplitudes
   * \f$T_{a}^{i}\f$ and \f$T_{ab}^{ij}\f$ and the Coulomb integrals \f$V_{ij}^{ab}\f$. For
   * calculating the amplitudes it calls the iteration routine of the actual algorithm.
   */
  class ClusterSinglesDoublesAlgorithm: public ClusterDoublesAlgorithm {
  public:
    ClusterSinglesDoublesAlgorithm(
      std::vector<Argument> const &argumentList
    );
    virtual ~ClusterSinglesDoublesAlgorithm();
    virtual void run();
    /**
     * \brief Returns the abbreviation of the concrete algorithm, e.g. "CCSD".
     */
    virtual std::string getAbbreviation() = 0;

  protected:
    /**
     * \brief The mixer for the singles amplitudes, additionally to those
     * of the inheritided doubles amplitudes TabijMixer
     */
    Mixer<double> *TaiMixer;

    /**
     * \brief Performs one iteration of the concrete algorithm.
     */
    virtual void iterate(int i) = 0;

    void singlesAmplitudesFromResiduum(CTF::Tensor<> &Rai);

    /**
     * \brief Calculates and returns one slice Xxycd of the Coulomb integrals
     * coupled to the singles amplitudes.
     * The indices x and y are restricted to the
     * range {No+a, ..., No+a+No-1} and {No+b, ..., No+b+No-1}, respectively.
     * The caller is responsible for deleting the dynamically allocated
     * result tensor. 
     */
    CTF::Tensor<> *sliceCoupledCoulombIntegrals(int a, int b, int sliceRank);
  };
}

#endif

