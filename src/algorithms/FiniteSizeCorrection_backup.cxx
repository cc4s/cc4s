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
  Tensor<complex> comVG(1, VG->lens, VG->sym, *VG->wrld, "comVG");
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
    Tensor<> *InvVG(getTensorArgument<>("invVG"));
    Tensor<complex> CGai(*GammaGai);
    Tensor<complex> comInvVG(1, VG->lens, VG->sym, *VG->wrld, "comInvVG");

    toComplexTensor(invVG, comInvVG);
    toComplexTensor(*VG, comVG);
    CGai["Gai"] *= comInvVG["G"];

    Tensor<complex> conjCGai(false, CGai);
    Univar_Function<complex> fConj(conj<complex>);
    conjCGai.sum(1.0,CGai,"Gai", 0.0,"Gai", fConj);

//    Tensor<> realCGai(3, CGai.lens, CGai.sym, *CGai.wrld, "realCGai");
//    Tensor<> imagCGai(3, CGai.lens, CGai.sym, *CGai.wrld, "imagCGai");
//    fromComplexTensor(CGai, realCGai, imagCGai);
//    toComplexTensor(realPart, imagPart, complexTensor);
//    toComplexTensor(realPart, complexTensor);

    Tensor<> *realTabij(getTensorArgument("DoublesAmplitudes"));
    Tensor<complex> Tabij(
      4, realTabij->lens, realTabij->sym, *realTabij->wrld, "Tabij"
    );
    toComplexTensor(*realTabij, Tabij);
    int NG(CGai.lens[0]);
    Vector<> *SG(new Vector<>(NG, *CGai.wrld, "SG"));
    Vector<complex> comSG(NG, *CGai.wrld, "comSG");
    // SG = Tabij ( 2*conj(CGai)*CGbj - conj(CGaj)*CGbi )
    comSG["G"] = 2.0 *  conjCGai["Gai"] * CGai["Gbj"] * Tabij["abij"];
    comSG["G"] -=  conjCGai["Gaj"] * CGai["Gbi"] * Tabij["abij"];
    Scalar<> EMp2(*CGai.wrld);
    fromComplexTensor(comSG, *SG);
    EMp2[""] = (*SG)["G"] * (*VG)["G"];
    setRealArgument("EMp2", EMp2);
    allocatedTensorArgument<>("StructureFactor", SG);
    allocatedTensorArgument<>("InvVG", InvVG);
   /* int NG(CGai.lens[0]);
    Vector<> *SG(new Vector<>(NG, *CGai.wrld, "SG"));
    (*SG)["G"] = 2.0 * realCGai["Gai"] * realCGai["Gbj"]* Tabij["abij"];
    (*SG)[""] +=2.0 * imagCGai["Gai"] * imagCGai["Gbj"]* Tabij["abij"];
    (*SG)["G"] -= realCGai["Gaj"] * realCGai["Gbi"]* Tabij["abij"];
    (*SG)["G"] -= imagCGai["Gaj"] * imagCGai["Gbi"]* Tabij["abij"];   
  */
   /*
    Tensor<complex> Tabij(
      4, realTabij->lens, realTabij->sym, *realTabij->wrld, "Tabij"
    );
    toComplexTensor(*realTabij, Tabij);
    Tensor<complex> Tempabij(false, Tabij);
    Tempabij["abij"] = 2.0 * conjCGai["Gai"] * CGai["Gbj"];
    Tempabij["abij"] -= conjCGai["Gaj"] * CGai["Gbi"];
   */
  }
  
  double energyCorrection(0.0);
  setRealArgument("EnergyCorrection", energyCorrection);
  allocatedTensorArgument<>("VG",VG);
}

