/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LEGACY_CLUSTER_SINGLES_DOUBLES_ALGORITHM_DEFINED 
#define LEGACY_CLUSTER_SINGLES_DOUBLES_ALGORITHM_DEFINED

#include <algorithms/LegacyClusterDoublesAlgorithm.hpp>
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
   */
  class LegacyClusterSinglesDoublesAlgorithm:
    public LegacyClusterDoublesAlgorithm
  {
  public:
    LegacyClusterSinglesDoublesAlgorithm(
      std::vector<Argument> const &argumentList
    );
    virtual ~LegacyClusterSinglesDoublesAlgorithm();
    /**
     * \brief Calculates the energy of a ClusterSinglesDoubles algorithm
     */
    virtual void run();
    /**
     * \brief Returns the abbreviation of the concrete algorithm in camel case,
     * e.g. "Ccsd", "Dcsd". They will constitute the first part of the
     * resulting tensors, such as "CcsdDoublesAmplitudes".
     */
    virtual std::string getAbbreviation() = 0;
    /**
     * \brief Performs a Dry Run
     */
    virtual void dryRun();

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

    /**
     * \brief Calculates the energy from the amplitudes currently contained
     * in the mixers. Overrides ClusterDoublesAlgorithm method.
     **/
    virtual double calculateEnergy();

    /**
     * \brief Calculates the singles amplitudes from the current residuum and
     * returns them in-place.
     * Usually this is done by calculating
     * \f$T_{i}^{a} = R_{i}^{a} / (\varepsilon_i-\varepsilon_a)\f$,
     * but other methods, such as level shifting may be used.
     * \param[in] Rai Residuum Tensor.
     */
    void singlesAmplitudesFromResiduum(CTF::Tensor<> &Rai);

    /**
     * \brief Dry run for singlesAmplitudesFromResiduum.
     * \param[in] Rai Residuum Tensor.
     */
    void drySinglesAmplitudesFromResiduum(cc4s::DryTensor<> &Rai);

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
    CTF::Tensor<> *sliceCoupledCoulombIntegrals(int a, int b, int integralsSliceSize);

    /**
     * \brief Dry run for sliceCoupledCoulombIntegrals. 
     * \param[in] a 1st sliced dimension (x).
     * \param[in] b 2nd sliced dimension (y).
     * \param[in] integralsSliceSize slicing rank.
     */
    cc4s::DryTensor<> *drySliceCoupledCoulombIntegrals(int integralsSliceSize);

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
    CTF::Tensor<> *sliceAmplitudesFromCoupledCoulombFactors(int a, int b, int factorsSliceSize);

    /**
     * \brief Calculates and returns tensor Fabij of the residuum
     * from the dressed Coulomb factors. 
     * The caller is responsible for deleting the dynamically allocated
     * result tensor. 
     * \param[out] Fabij sliced Residuum
     */
    CTF::Tensor<> *amplitudesFromCoupledCoulombFactors();

    /**
     * \brief Dry run for sliceAmplitudesFromCoupledCoulombFactors.
     * \param[in] factorsSliceSize slicing rank of NR.
     * \param[out] Fabij sliced Residuum
     */
    cc4s::DryTensor<> *drySliceAmplitudesFromCoupledCoulombFactors(int factorsSliceSize);
  };
}

#endif

