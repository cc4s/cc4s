/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COULOMB_VERTEX_SINGULAR_VECTORS_DEFINED
#define COULOMB_VERTEX_SINGULAR_VECTORS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class CoulombVertexSingularVectors: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CoulombVertexSingularVectors);
    CoulombVertexSingularVectors(std::vector<Argument> const &argumentList);
    virtual ~CoulombVertexSingularVectors();
    /**
     * \brief Calculates left singular vectors of the full Coulomb vertex
     * \f$\tilde\Gamma^q_{rG}\f$.
     */
    virtual void run();
    /**
     * \brief Dry run for calculating the left singular vectors of
     * \f$\Gamma^q_{rG}\f$.
     */
    virtual void dryRun();

    static double constexpr DEFAULT_FIELD_VARIABLES_RANK = 0.5;
    static int64_t constexpr DEFAULT_FIELD_VARIABLES_SIZE = -1;
  };
}

#endif

