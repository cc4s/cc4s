/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DOUBLES_AMPLITUDES_FROM_VERTEX_DEFINED
#define DOUBLES_AMPLITUDES_FROM_VERTEX_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  /**
   * \brief Caclulates the doubles amplitudes \f$T_{ij}^{ab}\f$
   * from the doubles amplitudes Vertex \f$Y^a_{iL}\f$.
   */
  class DoublesAmplitudesFromVertex: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DoublesAmplitudesFromVertex);
    DoublesAmplitudesFromVertex(
      std::vector<Argument> const &argumentList
    );
    virtual ~DoublesAmplitudesFromVertex();
    /**
     * \brief Calculates doubles amplitudes Tabij.
    */
    virtual void run();
    /**
     * \brief Dry run for calculating doubles integrals Tabij
    */
    virtual void dryRun();

  };
}

#endif

