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

void ThermalDirectRingCoupledClusterDoubles::applyHamiltonian(
  const CTF::Tensor<real> &T0abij,
  CTF::Tensor<real> &T1abij,
  real DTau
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

