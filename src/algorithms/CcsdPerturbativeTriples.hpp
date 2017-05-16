/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CCSD_PERTURBATIVE_TRIPLES_DEFINED
#define CCSD_PERTURBATIVE_TRIPLES_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/Permutation.hpp>
#include <util/SlicedCtfTensor.hpp>

namespace cc4s {
  /**
   * \brief Caclulates perturbative triples correction
   */
  class CcsdPerturbativeTriples: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CcsdPerturbativeTriples);
    CcsdPerturbativeTriples(std::vector<Argument> const &argumentList);
    virtual ~CcsdPerturbativeTriples();
    /**
     * \brief Calculates perturbative triples correction. Routine based on Helgaker book.
     */
    virtual void run();

    /**
     * \brief Dry run for perturbative triples correction based on Helgaker book.
     */
    virtual void dryRun();

  protected:
    int No, Nv;
    CTF::Tensor<> *SVabc, *DVabc;
    CTF::Tensor<> *realGammaFab, *imagGammaFab;
    SlicedCtfTensor<> *Tai, *Tabij, *Tabil;
    SlicedCtfTensor<> *Vabij, *Vijla, *realGammaFai,*imagGammaFai;
    void sliceTensors();
    CTF::Tensor<> &getSinglesContribution(const Map<3> &);
    CTF::Tensor<> &getDoublesContribution(const Map<3> &);
    CTF::Tensor<> &getEnergyDenominator(const Map<3> &);
  };
}

#endif

