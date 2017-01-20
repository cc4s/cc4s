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
     * \brief Calculates perturbative triples correction. Routine based on Helgaker book.
     */
    virtual void runHelgaker();
    /**
     * \brief Calculates perturbative triples correction. Routine based on Piecuch paper.
     */
    virtual void runPiecuch();
    /**
     * \brief Calculates perturbative triples correction. Routine based on Piecuch paper.
     */
    virtual void run();
    /**
     * \brief Dry run for perturbative triples correction based on Helgaker book.
     */
    virtual void dryRunHelgaker();
    /**
     * \brief Dry run for perturbative triples correction based on Piecuch paper.
     */
    virtual void dryRun();
  };
}

#endif

