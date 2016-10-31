/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PERTURBATIVE_TRIPLES_DEFINED
#define PERTURBATIVE_TRIPLES_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  /**
   * \brief Caclulates perturbative triples correction
   */
  class PerturbativeTriples: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(PerturbativeTriples);
    PerturbativeTriples(std::vector<Argument> const &argumentList);
    virtual ~PerturbativeTriples();
    /**
     * \brief Calculates perturbative triples correction
     */
    virtual void run();
    /**
     * \brief Dry run for perturbative triples correction
     */
    virtual void dryRun();
  };
}

#endif

