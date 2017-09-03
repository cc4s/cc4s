#include <algorithms/ThermalDirectRingCoupledClusterDoubles.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ThermalDirectRingCoupledClusterDoubles);

ThermalDirectRingCoupledClusterDoubles::ThermalDirectRingCoupledClusterDoubles(
  std::vector<Argument> const &argumentList
): ThermalClusterDoublesAlgorithm(argumentList) {
}

ThermalDirectRingCoupledClusterDoubles::
  ~ThermalDirectRingCoupledClusterDoubles(
) {
}

void ThermalDirectRingCoupledClusterDoubles::iterate(int n) {
  // Read Vabij
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));

  Tensor<> PVabij(*Vabij);

  // get the occupancies for contractions
  Tensor<> *Ni(getTensorArgument("ThermalHoleOccupancies"));
  Tensor<> *Na(getTensorArgument("ThermalParticleOccupancies"));

  PVabij["abij"] *= (*Ni)["i"];
  PVabij["abij"] *= (*Ni)["j"];
  PVabij["abij"] *= (*Na)["a"];
  PVabij["abij"] *= (*Na)["b"];
}


void ThermalDirectRingCoupledClusterDoubles::dryIterate() {
  // Read the DRCCD amplitudes Tabij
  //DryTensor<> *Tabij(
  getTensorArgument<double, DryTensor<double>>("DrccdDoublesAmplitudes");
  //);

  // Read the Coulomb Integrals Vabij
  DryTensor<> *Vabij(getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals"));

  // Allocate Tensors for T2 amplitudes
  DryTensor<> Rabij(*Vabij);

  // Define intermediates
  DryTensor<> Cabij(*Vabij);

  DryTensor<> Dabij(*Vabij);
}

