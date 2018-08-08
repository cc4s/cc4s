/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef FINITE_SIZE_CORRECTION_DEFINED
#define FINITE_SIZE_CORRECTION_DEFINED

#include <algorithms/Algorithm.hpp>
#include <math/Interpolation.hpp>
#include <math/Vector.hpp>
namespace cc4s {
  class FiniteSizeCorrection: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(FiniteSizeCorrection);
    FiniteSizeCorrection(
      std::vector<Argument> const &argumentList
    );
    virtual ~FiniteSizeCorrection();

    /**
     * \brief Calculates the finite size correction
     */
    virtual void run();

    /**
     * \brief Performs a Dry Run
     */
    virtual void dryRun();
  protected:
    int NG;
    static constexpr int DEFAULT_NUM_FIBONACCI = 128;
    std::vector<double> GLengths;
    std::vector<double> averageSGs;
    std::vector<double> meanErrorSG;
    std::vector<double> structureFactors;
    std::vector<double> VofG;
    double GC;
    double inter3D;
    double sum3D;
    class Momentum;
    Momentum *fibonacciGrid;
    Momentum *cartesianGrid;
    void readFromFile();
    void calculateRealStructureFactor();
    void calculateRealStructureFactorSliced(int iStart, int iEnd);
    //    template <typename F>
    //    F calculateStructureFactor(CTF::Tensor<F> &Vabij);

    void calculateComplexStructureFactor();

    void constructFibonacciGrid(double R, int N);
    void interpolation3D();
    bool IsInSmallBZ(
      Vector<double> point, double scale, std::vector<cc4s::Vector<double>> smallBZ
    );
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

    void dryCalculateStructureFactor();
    void dryInterpolation3D();
    void dryCalculateFiniteSizeCorrection();

    void extrapolation(double minG, double maxG, int basisSetExtrapolation);
    double simplestWindow(double Gmin, double Gmax, double G);
    double integrateSimplestWindow(double Gmin, double Gmax);
    double leastSquareFit(std::vector<double> fitabsG, std::vector<double> fitSF);

    void basisSetCompleteness();
 
 };
}

#endif
