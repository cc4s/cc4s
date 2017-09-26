/*Copyright (c) 2017, Andreas Grueneis, Felix Hummel and Alejandro Gallo, all rights reserved.*/
#ifndef MP2_EOM_DEFINED
#define MP2_EOM_DEFINED

#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>

namespace cc4s {
  /**
   * \brief Implements the iteration routine for the Mp2 method. Calculates the
   * amplitudes \f$T_{ab}^{ij}\f$ from the Coulomb Integrals \f$V_{ij}^{ab}\f$
   * in a \f$ \mathcal{O}(N^{6}) \f$ implementation.
   */
  class UccsdAmplitudesFromCoulombIntegrals: public ClusterSinglesDoublesAlgorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(UccsdAmplitudesFromCoulombIntegrals);
    UccsdAmplitudesFromCoulombIntegrals(
      std::vector<Argument> const &argumentList
    );
    virtual ~UccsdAmplitudesFromCoulombIntegrals();

    virtual void run();
    double calculateEnergy();
    virtual std::string getAbbreviation() { return "Uccsd"; }

  protected:
    virtual void iterate(int i);
    CTF::Tensor<> *epsi;
    CTF::Tensor<> *epsa;

    // Get couloumb integrals (these shoul not be antisymmetrized)
    CTF::Tensor<> *Vijkl;
    CTF::Tensor<> *Vabcd;
    CTF::Tensor<> *Vabij;
    CTF::Tensor<> *Vijka;
    CTF::Tensor<> *Vaibj;
    CTF::Tensor<> *Vabci;

    //  Vijab
    CTF::Tensor<> *Vijab;
    //  Viajk
    CTF::Tensor<> *Viajk;
    // Viajb
    CTF::Tensor<> *Viajb;
    // Viabc
    CTF::Tensor<> *Viabc;
    // Vabic
    CTF::Tensor<> *Vabic;



  };
}

#endif

