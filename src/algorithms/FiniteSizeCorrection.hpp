/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef FINITE_SIZE_CORRECTION_DEFINED
#define FINITE_SIZE_CORRECTION_DEFINED

#include <algorithms/Algorithm.hpp>

namespace cc4s {
  class FiniteSizeCorrection: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(FiniteSizeCorrection);
    FiniteSizeCorrection(
      std::vector<Argument> const &argumentList
    );
    virtual ~FiniteSizeCorrection();
    virtual void run();
  protected:
    int NG;
    double *structureFactors;

    void calculateStructureFactor();
    void constructFibonacciGrid();
    void interpolation3D();
    void calculateFiniteSizeCorrection();
  };
}

#endif

