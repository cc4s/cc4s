/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef DOUBLES_AMPLITUDES_SINGULAR_VECTORS_DEFINED
#define DOUBLES_AMPLITUDES_SINGULAR_VECTORS_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class DoublesAmplitudesSingularVectors: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(DoublesAmplitudesSingularVectors);
    DoublesAmplitudesSingularVectors(std::vector<Argument> const &argumentList);
    virtual ~DoublesAmplitudesSingularVectors();
    /**
     * \brief Calculates left singular vectors of the particle-hole Coulomb vertex
     * \f$\tilde\Gamma^a_{iG}\f$.
     */
    virtual void run();
    /**
     * \brief Dry run for calculating the left singular vectors of
     * \f$\Gamma^a_{iG}\f$.
     */
    virtual void dryRun();

    static double constexpr DEFAULT_REDUCTION = 2.0;
    static int64_t constexpr DEFAULT_FIELD_VARIABLES = -1;
  };
}

#endif

