#include <algorithms/FiniteSizeCorrection.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <util/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(FiniteSizeCorrection);

FiniteSizeCorrection::FiniteSizeCorrection(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

FiniteSizeCorrection::~FiniteSizeCorrection() {
}

void FiniteSizeCorrection::run() {
  Tensor<complex> *GammaGai(
    getTensorArgument<complex>("ParticleHoleCoulombVertex")
  );
  Tensor<> *VG(getTensorArgument<>("CoulombKernel"));
  Tensor<> invVG(false, *VG);
  {
    class InvSqrt {
    public:
      double operator ()(double x) {
//        LOG_RANK(1, "FiniteSizeCorrection") << "is inf? " << std::isinf(x) <<"invert of "<< std::sqrt(1.0/x) << std::endl;
        return std::sqrt(1.0 / x);
      }
    };
    InvSqrt invSqrt;
    Univar_Function<> fInvSqrt(invSqrt);
    invVG.sum(1.0,*VG,"G", 0.0,"G", fInvSqrt);
    Tensor<complex> CGai(*GammaGai);
    CGai["Gai"] *= invVG["G"];

    Tensor<complex> conjCGai(false, CGai);
    Univar_Function<complex> fConj(conj<complex>);
    conjCGai.sum(1.0,CGai,"Gai", 0.0,"Gai", fConj);
/*
    Tensor<> realCGai(3, CGai.lens, CGai.sym, *CGai.wrld, "realCGai");
    Tensor<> imagCGai(3, CGai.lens, CGai.sym, *CGai.wrld, "realCGai");
    fromComplexTensor(CGai, realCGai, imagCGai);
*/
    Tensor<> *Tabij(getTensorArgument("DoublesAmplitudes"));
    int NG(CGai.lens[0]);
    Vector<> SG(NG, *CGai.wrld, "SG");
    // SG = Tabij ( 2*conj(CGai)*CGbj - conj(CGaj)*CGbi )
    SG["G"] = conjCGai["Gai"] * CGai["Gbj"] * (*Tabij)["abij"];
  }
  
  double energyCorrection(0.0);
  setRealArgument("EnergyCorrection", energyCorrection);
}

