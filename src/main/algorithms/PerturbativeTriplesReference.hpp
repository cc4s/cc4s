/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PERTURBATIVE_TRIPLES_DEFINED
#define PERTURBATIVE_TRIPLES_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  /**
   * \brief Reference implementation of the perturbative triples energy
   */
  class PerturbativeTriplesReference: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(PerturbativeTriplesReference)
    /**
     * \brief run routine
     */
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
  protected:
    template <typename F, typename TE>
    Ptr<MapNode> calculateTriplesEnergy(const Ptr<MapNode> &arguments);

  };
}

#endif

