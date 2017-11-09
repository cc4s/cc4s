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

  private:
    // NOTE: the Dummy template argument is needed to "fully" specialize
    // the inner class CoulombVertex to CoulombVertex<double> or <complex>
    template <typename F, int Dummy=0>
    class CoulombVertex {
    };

    template <int Dummy>
    class CoulombVertex<double, Dummy> {
    public:
      CoulombVertex(
        CTF::Tensor<complex> *GammaFab,
        CTF::Tensor<complex> *GammaFai
      );
      ~CoulombVertex();
      void getDoublesParticleContribution(
        SlicedCtfTensor<double> &Tabij, const Map<3> &i,
        CTF::Tensor<double> &SVabc
      );
    protected:
      CTF::Tensor<double> *realGammaFab, *imagGammaFab;
      SlicedCtfTensor<double> *realGammaFai,*imagGammaFai;
    };

    template <int Dummy>
    class CoulombVertex<complex, Dummy> {
    public:
      CoulombVertex(
        CTF::Tensor<complex> *GammaFab,
        CTF::Tensor<complex> *GammaFai
      );
      ~CoulombVertex();
      void getDoublesParticleContribution(
        SlicedCtfTensor<complex> &Tabij, const Map<3> &i,
        CTF::Tensor<complex> &DVabc
      );
    protected:
      CTF::Tensor<complex> *conjGammaFab;
      SlicedCtfTensor<complex> *GammaFai;
    };

    template <typename F>
    class Calculator {
    public:
      Calculator(
        CTF::Tensor<F> *Tai, CTF::Tensor<F> *Tabij,
        CTF::Tensor<F> *Vabij, CTF::Tensor<F> *Valij,
        CTF::Tensor<complex> *GammaFab,
        CTF::Tensor<complex> *GammaFai,
        CTF::Tensor<double> *epsi, CTF::Tensor<double> *epsa
      );
      ~Calculator();
      F calculate();
      void addDoublesHoleContribution(const Map<3> &, CTF::Tensor<F> &);
      CTF::Tensor<F> &getSinglesContribution(const Map<3> &);
      CTF::Tensor<F> &getEnergyDenominator(const Map<3> &);
    protected:
      CTF::Tensor<F> *SVabc, *DVabc;
      SlicedCtfTensor<F> *Tai, *Tabij, *Tabil;
      SlicedCtfTensor<F> *Vabij, *Valij;
      CoulombVertex<F> Gamma;
      SlicedCtfTensor<F> *epsi;
      CTF::Tensor<F> *epsa;
    };
  };
}

#endif

