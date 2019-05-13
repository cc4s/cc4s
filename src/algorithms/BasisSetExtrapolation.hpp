/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef BASIS_SET_EXTRAPOLATION_FUNCTION_DEFINED
#define BASIS_SET_EXTRAPOLATION_FUNCTION_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/Vector.hpp>
#include <vector>

namespace cc4s {
  class BasisSetExtrapolation: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(BasisSetExtrapolation);
    BasisSetExtrapolation(
      std::vector<Argument> const &argumentList
    );
    virtual ~BasisSetExtrapolation();
    virtual void run();
  protected:

    static constexpr int iterations = 10;

    void evaluateQGG(int slice);
    void fitF12(int type, real minG, real maxG);
    void calculateNewSF(
      int type, real gamma, CTF::Tensor<> *coulombKernel, CTF::Tensor<> *newSF, CTF::Tensor<> *resNewSF
    );
    void invertQGG();
  };
}


#endif
