/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef REDUCE_COULOMB_VERTEX_DEFINED
#define REDUCE_COULOMB_VERTEX_DEFINED

#include <algorithms/Algorithm.hpp>
#include <ctf.hpp>

namespace cc4s {
  class ReduceCoulombVertex: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(ReduceCoulombVertex);
    ReduceCoulombVertex(std::vector<Argument> const &argumentList);
    virtual ~ReduceCoulombVertex();
    /**
     * \brief Reduces the Coulomb vertex by applying the
     * energy matrix reduction transform:
     * \f$\Gamma^{qg}_r = \Gamma^{qG}_r U_G^g\f$.
     */
    virtual void run();
    /**
     * \brief Dry run of reducing the Coulomb vertex.
     */
    virtual void dryRun();
  };
}

#endif

