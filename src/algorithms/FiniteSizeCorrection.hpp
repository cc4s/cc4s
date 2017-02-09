/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef FINITE_SIZE_CORRECTION_DEFINED
#define FINITE_SIZE_CORRECTION_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/Interpolation.hpp>
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
    const int  N=128;//The # of points on the fibonacciGrid, fixed #.
    int numBins;
    double *lengthG;
    double *avgSG;
    double *structureFactors;
    double *VofG;
    double GC;
    class Momentum;
    Momentum *fibonacciGrid;
    Momentum *momentumGrid;
    void calculateStructureFactor();
    void constructFibonacciGrid(double R);
    void interpolation3D();
    double SGxVG(cc4s::Inter1D<double> Int1d, double x);
    double integrate(
      cc4s::Inter1D<double> Int1d,
      double start, double end, int steps
      );
    double simpson(
      cc4s::Inter1D<double> Int1d,
      double x, double h
      );
    void calculateFiniteSizeCorrection();
  };
}

#endif

