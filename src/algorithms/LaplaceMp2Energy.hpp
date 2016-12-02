/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LAPLACE_MP2_ENERGY_DEFINED
#define LAPLACE_MP2_ENERGY_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/Complex.hpp>

namespace cc4s {
  /**
   * \brief Caclulates MP2 energy from the Coulomb Integrals \f$V_{ij}^{ab}.
   */
  class LaplaceMp2Energy: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(LaplaceMp2Energy);
    LaplaceMp2Energy(std::vector<Argument> const &argumentList);
    virtual ~LaplaceMp2Energy();
    /**
     * \brief Calculates MP2 energy from Coulomb integrals Vabij
     */
    virtual void run();
    /**
     * \brief Dry run for the MP2 energy from Coulomb integrals Vabij
     */
    virtual void dryRun();

  protected:
    CTF::Tensor<complex> *GpRSn, *GhRSn, *VRS, *wn;

    double calculateDirectTerm();
    double calculateExchangeTerm();
    double calculateDirectTermAnalytically();
  };
}

#endif

