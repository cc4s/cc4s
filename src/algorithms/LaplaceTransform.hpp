/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LAPLACE_TRANSFORM_DEFINED 
#define LAPLACE_TRANSFORM_DEFINED

#include <algorithms/Algorithm.hpp>
#include <util/SharedPointer.hpp>
#include <string>
#include <ctf.hpp>
#include <tcc/DryTensor.hpp>

namespace cc4s {
  /**
   * \brief Tests the numerical Laplace transformation grids for the
   * propagators exp(-Delta*tau), where Delta=eps_a-eps_i.
   * Note that in a thermal system Delta can be negative.
   **/
  class LaplaceTransform: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(LaplaceTransform);
    LaplaceTransform(
      std::vector<Argument> const &argumentList
    );
    virtual ~LaplaceTransform();
    /**
     * \brief Writes the forward and backward transformation errors.
     */
    virtual void run();

  protected:
    std::string getAmplitudeIndices(CTF::Tensor<> &T);
    void fetchDelta(CTF::Tensor<> &Delta);
  };
}

#endif

