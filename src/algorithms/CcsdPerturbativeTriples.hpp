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
    template <typename F>
    class Calculator;

    template <>
    class Calculator<double> {
    public:
      int No, Nv;
      
      CTF::Tensor<> *realGammaFab, *imagGammaFab;
      SlicedCtfTensor<> *realGammaFai,*imagGammaFai;

      //      CTF::Tensor<complex> *GammaFab;
      //      SlicedCtfTensor<complex> *GammaFai;

      CTF::Tensor<F> *SVabc, *DVabc;
      SlicedCtfTensor<F> *Tai, *Tabij, *Tabil;
      SlicedCtfTensor<F> *Vabij, *Vijla;

      F run();

      void sliceTensors();
      CTF::Tensor<F> &getSinglesContribution(const Map<3> &);
      CTF::Tensor<F> &getDoublesContribution(const Map<3> &);
      CTF::Tensor<F> &getEnergyDenominator(const Map<3> &);
    };
  };
}

#endif

