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
    virtual void dryRun();
  protected:
    int NG;
    int NR;
    int No;
    int Nv;
    int Np;
    int NF;
    std::vector<real> structureFactor;
    std::vector<real> pairCorrelationFunction;
    std::vector<Vector<>> reciprocalGrid;
    std::vector<Vector<>> realSpaceMesh;
    std::vector<real> realSpaceBasisSetCompleteness;    
    
    void extrapolation(real Gmin, real Gmax, int itype);
    void basisSetCompletenessFullFTODDUMP();
    void basisSetCompletenessFTODDUMPIA();
    real leastSquareFit(std::vector<real> fitabsG, std::vector<real> fitSF);
    real simplestWindow(real Gmin, real Gmax, real G);
    real integrateSimplestWindow( real Gmin, real Gmax);
    void GtoRFourier();
    void RtoGFourier();
    void readReciprocalGridFromFile();
    void constructRealSpaceMesh(int augfac,std::vector<Vector<>>& realSpaceMesh);
    void fourierCompleteness();
    void inverseConvolution();
    void QGG(int iStart, int iEnd);
    void dryQGG(int iStart, int iEnd);
    void realCGi();
    void fitF12(int type, int iter, real minG, real maxG);
    void calculateNewSF(int type, real gamma, CTF::Tensor<> *coulombKernel, CTF::Tensor<> *newSF, CTF::Tensor<> *resNewSF);
  };
}


#endif
