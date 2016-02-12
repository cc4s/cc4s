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
  Tensor<> *Vabij(getTensorArgument("ParticleHoleCoulombIntegrals"));
 
  Tensor<> Tabij(false, Vabij);
  Tabij["abij"] =  (*epsi)["i"];
  Tabij["abij"] += (*epsi)["j"];
  Tabij["abij"] -= (*epsa)["a"];
  Tabij["abij"] -= (*epsa)["b"];

  Bivar_Function<> fDivide(&divide<double>);
  Tabij.contract(1.0, (*Vabij),"abij", Tabij,"abij", 0.0,"abij", fDivide);

  Scalar<> energy(*Cc4s::world);
  double e, dire, exce;

  energy[""] = 2.0 * Tabij["abij"] * (*Vabij)["abij"];
  dire = energy.get_val();
  energy[""] = Tabij["abji"] * (*Vabij)["abij"];
  exce = -1.0 * energy.get_val();
  e = dire + exce;

  LOG(0) << "e=" << e << std::endl;
  LOG(1) << "MP2d=" << dire << std::endl;
  LOG(1) << "MP2x=" << exce << std::endl;

  setRealArgument("Mp2Energy", e);
}

