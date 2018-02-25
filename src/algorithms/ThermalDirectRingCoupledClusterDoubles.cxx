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

void ThermalDirectRingCoupledClusterDoubles::getResiduum(
  CTF::Tensor<complex> &Tabijn
) {
/*
  // * interaction V
  Tabij["acik"] *= (*Vabij)["acik"];
  // * thermal weight of contracted indices ck
  Fck["ck"] = 1.0;
  thermalContraction(Fck);
  Wabij["acik"] *= Fck["ck"];
*/
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

