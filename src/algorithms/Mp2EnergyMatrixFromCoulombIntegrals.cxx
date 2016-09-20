#include <algorithms/Mp2EnergyMatrixFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <util/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <cmath>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(Mp2EnergyMatrixFromCoulombIntegrals);

Mp2EnergyMatrixFromCoulombIntegrals::Mp2EnergyMatrixFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

Mp2EnergyMatrixFromCoulombIntegrals::~Mp2EnergyMatrixFromCoulombIntegrals() {
}

/*
class raise {
public:
  raise(double expon_): expon(expon_) {
  }
  double operator()(double element) {
    //    return std::pow(element,expon);
    return 0;
  }
protected:
  double expon;
};
*/

void Mp2EnergyMatrixFromCoulombIntegrals::run() {
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
 
  Tensor<> Tabij(Vabij);
  Tabij["abij"] *= 2.0;
  Tabij["abij"] -= (*Vabij)["abji"];

  Tensor<> Dabij(false, *Vabij);
  Dabij["abij"] -= (*epsi)["i"];
  Dabij["abij"] -= (*epsi)["j"];
  Dabij["abij"] += (*epsa)["a"];
  Dabij["abij"] += (*epsa)["b"];

  //  double exponent(getRealArgument("exponent" , 1.0));
    
  //  Univar_Function<> raiseTensor(raise(exponent));
  //  Univar_Function<> raiseTensor(&sqrt<double>);
  //  Tensor<> raiseDabij(false, Dabij);
  //  raiseDabij.sum(1.0, Dabij, "abij", 0.0, "abij", &raiseTensor);
    
  Bivar_Function<> fDivide(&divide<double>);
  Tabij.contract(-1.0, Tabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);

  // TODO: use complex conversion routines
  Tensor<> imagTabij(false, Tabij);
  Tensor<complex> cTabij(4, Tabij.lens, Tabij.sym, *Tabij.wrld, "cTabij");
  // convert into complex tensor
  toComplexTensor(Tabij, imagTabij, cTabij);

  Tensor<complex> *GammaGai(getTensorArgument<complex>("ParticleHoleCoulombVertex"));
  // calculate conjugate of GammaGai
  Tensor<complex> conjGammaGai(*GammaGai);
  conjugate(conjGammaGai);

  Matrix<complex> *energyMatrix = new Matrix<complex>(
    GammaGai->lens[0], GammaGai->lens[0], *Vabij->wrld
  );
  allocatedTensorArgument<complex>("Mp2EnergyMatrix", energyMatrix);

  (*energyMatrix)["GH"] = conjGammaGai["Gai"] * (*GammaGai)["Hbj"] * cTabij["abij"];
}

void Mp2EnergyMatrixFromCoulombIntegrals::dryRun() {
  //DryTensor<> *Vabij(
  getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals");
  //);

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies")
  );
  
  // Compute the No,Nv
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Allocate the doubles amplitudes
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { Nv, Nv, No, No };
  DryTensor<> Tabij(4, vvoo, syms);

  DryTensor<complex> *GammaGai(
    getTensorArgument<complex, DryTensor<complex>>("ParticleHoleCoulombVertex")
  );

  DryMatrix<complex> *energyMatrix = new DryMatrix<complex>(
    GammaGai->lens[0], GammaGai->lens[0], NS
  );
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "Mp2EnergyMatrix", energyMatrix
  );
}

