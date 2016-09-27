#include <algorithms/FiniteSizeCorrection.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
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
  {
    Tensor<complex> CGai(GammaGai);
    class InvSqrt {
    public:
      double operator ()(double x) {
        LOG_RANK(1, "FiniteSizeCorrection") << "invSqrt of " << x << std::endl;
        return std::isinf(x) ? 0.0 : std::sqrt(1.0 / x);
      }
    };
    InvSqrt invSqrt;
    Univar_Function<> fInvSqrt(invSqrt);
    VG->sum(1.0,*VG,"G", 0.0,"G", fInvSqrt);
  }
  
  double energyCorrection(0.0);
  setRealArgument("EnergyCorrection", energyCorrection);
}

