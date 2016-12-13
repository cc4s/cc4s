/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COULOMB_VERTEX_FROM_FACTORS_DEFINED
#define COULOMB_VERTEX_FROM_FACTORS_DEFINED

#include <algorithms/Algorithm.hpp>
#include <tcc/MachineTensor.hpp>
#include <memory>

namespace cc4s {
  /**
   * \brief Caclulates the Coulomb vertex \f$\Gamma^q_{rF}\f$ from the given
   * given factor orbitals \f$\Pi^R_r\f$ and Coulomb factors \f$\Lambda^R_F\f$.
   */
  class CoulombVertexFromFactors: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CoulombVertexFromFactors);
    CoulombVertexFromFactors(
      std::vector<Argument> const &argumentList
    );
    virtual ~CoulombVertexFromFactors();

    virtual void run();
    virtual void dryRun();

    template <typename T, typename MT>
    void runIt(const bool dryRun);
  };
}

#endif

