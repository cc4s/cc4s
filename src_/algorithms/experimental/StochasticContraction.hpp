/*Copyright (c) 2017, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef STOCHASTIC_CONTRACTION_DEFINED
#define STOCHASTIC_CONTRACTION_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/RandomGenerator.hpp>

namespace cc4s {
  /**
   * \brief Determines variances of distributions numerically for
   * stocastically evaluating tensor contractions
   */
  class StochasticContraction: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(StochasticContraction);
    StochasticContraction(std::vector<Argument> const &argumentList);
    virtual ~StochasticContraction();
    virtual void run();
  private:
    RandomGenerator *rand;
    complex drawUniformWeight();
  };
}

#endif

