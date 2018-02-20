/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef RPA_APX_ENERGY_DEFINED 
#define RPA_APX_ENERGY_DEFINED

#include <algorithms/Algorithm.hpp>

#include <util/SharedPointer.hpp>
#include <tcc/Tcc.hpp>

namespace cc4s {
  class RpaApxEnergy: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(RpaApxEnergy);
    RpaApxEnergy(
      std::vector<Argument> const &argumentList
    );
    virtual ~RpaApxEnergy();

    virtual void run();
  protected:
    void diagonalizeChiV();

    /**
     * \brief \f$(V^{1/2}.\chi.V^{1/2})^F_G\f$ for each
     * imaginary freqency point \f$\nu_n\f$.
     * After diagonalization, it contains its eigenvectors \f$U^L_G\f$ such
     * than \f$(V^{1/2}.\chi.V^{1/2})^F_G = {U^\ast}^F_L (\chi V)_L U^L_G\f$
     * where \f$(\chi V)_L\f$ the eigenvalues of
     * \f$(V^{1/2}.\chi.V^{1/2})^F_G\f$.
     **/
    PTR(tcc::Tensor<complex>) chi0VFGn, chi1VFGn;
    /**
     * \brief \f$(\chi V)_L\f$ are the eigenvalues of
     * \f$(V^{1/2}.\chi.V^{1/2})^F_G\f$ for each imaginary frequency point
     * \f%\nu_n\f$.
     **/
    PTR(tcc::Tensor<complex>) chi0VLn;
  };
}

#endif

