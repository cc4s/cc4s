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
//Definition of the variables
  Tensor<complex> *GammaGai(
	getTensorArgument<complex>("ParticleHoleCoulombVertex")
  );
  Tensor<> *VG(getTensorArgument<>("CoulombKernel"));
  Tensor<> realInvSqrtVG(false, *VG);
  Tensor<complex> invSqrtVG(
    1, realInvSqrtVG.lens, realInvSqrtVG.sym, 
     *realInvSqrtVG.wrld, "invSqrtVG"
  );
//Starting a new space whose memory will be erased after operation
  
    //Define operation inverse square root
    class InvSqrt {
    public:
	double operator ()(double x){
		return std::sqrt(1.0 / x);
	}
    };
  //Get the inverted square root of VG
  InvSqrt invSqrt;
  Univar_Function<> fInvSqrt(invSqrt);
  realInvSqrtVG.sum(1.0, *VG, "G", 0.0, "G", fInvSqrt);
  toComplexTensor(realInvSqrtVG, invSqrtVG);

  //Define CGai
  Tensor<complex> CGai(*GammaGai);
//  CGai["Gai"] *= invSqrtVG["G"];

  //Conjugate of CGai
  Tensor<complex> conjCGai(false, CGai);
  Univar_Function<complex> fConj(conj<complex>);
  conjCGai.sum(1.0, CGai, "Gai", 0.0, "Gai", fConj);

  //Get Tabij
  Tensor<> *realTabij(getTensorArgument("DoublesAmplitudes"));
  Tensor<complex> Tabij(
    4, realTabij->lens, realTabij->sym, *realTabij->wrld, "Tabij"
  );
  toComplexTensor(*realTabij, Tabij);
  
  //construct SG
  int NG(CGai.lens[0]);
  Vector<complex> *SG(new Vector<complex>(NG, *CGai.wrld, "SG"));
  (*SG)["G"] = 2.0 * conjCGai["Gai"] * CGai["Gbj"] * Tabij["abij"];
  (*SG)["G"] -= conjCGai["Gaj"] * CGai["Gbi"] * Tabij["abij"];
  Vector<> *realSG(new Vector<>(NG, *CGai.wrld, "realSG"));
  fromComplexTensor(*SG, *realSG);
  allocatedTensorArgument<>("StructureFactor", realSG);
  //Get EMp2
  Scalar<complex> EMp2(*CGai.wrld);
  EMp2[""] = (*SG)["G"]; // * (*VG)["G"];
  setRealArgument("EMp2", std::real(EMp2.get_val()));  
  
  
  allocatedTensorArgument<>("VG", VG);

}










