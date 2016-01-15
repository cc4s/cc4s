#include <Mp2EnergyFromCoulombIntegrals.hpp>
#include <util/Log.hpp>
#include <util/MathFunctions.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(Mp2EnergyFromCoulombIntegrals);

Mp2EnergyFromCoulombIntegrals::Mp2EnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

Mp2EnergyFromCoulombIntegrals::~Mp2EnergyFromCoulombIntegrals() {
}

/**
 * \brief Calculates MP2 energy from Coulomb integrals Vabij
 */
void Mp2EnergyFromCoulombIntegrals::run() {
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<> *vabij(getTensorArgument("ParticleHoleCoulombIntegrals"));
 
  Tensor<> Dabij(false, vabij);
  Tensor<> Tabij(false, vabij);

  Dabij["abij"] =  (*epsi)["i"];
  Dabij["abij"] += (*epsi)["j"];
  Dabij["abij"] -= (*epsa)["a"];
  Dabij["abij"] -= (*epsa)["b"];

  Bivar_Function<> fDivide(&divide<double>);
  Tabij.contract(1.0, (*vabij),"abij", Dabij,"abij", 0.0,"abij", fDivide);

  Scalar<> energy(*Cc4s::world);
  double e, dire, exce;

  energy[""] = 2.0 * Tabij["abij"] * (*vabij)["abij"];
  dire = energy.get_val();
  energy[""] = Tabij["abji"] * (*vabij)["abij"];
  exce = -1.0 * energy.get_val();
  e = dire + exce;
  LOG(0) << "e=" << e << std::endl;
  LOG(1) << "MP2d=" << dire << std::endl;
  LOG(1) << "MP2x=" << exce << std::endl;

  setRealArgument("Mp2Energy", e);
}

