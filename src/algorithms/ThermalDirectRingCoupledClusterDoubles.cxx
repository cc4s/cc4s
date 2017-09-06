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
  const int m(n+recursionLength);
  const double betam( beta*std::pow(recursionScaling,-m) );

  // Read Vabij
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  Tensor<> nextT( *amplitudes[recursionLength] );
  Tensor<> Tabij(false, *Vabij);
  Tensor<> Wabij(false, *Vabij);
  Tensor<> Dai(2, &Vabij->lens[1], &Vabij->sym[1], *Vabij->wrld, "Dai");
  Tensor<> Fck(false, Dai);

  // V*T:
  // effective V = W = propagation of states connected to V
  fetchDelta(Dai);
  // states ck are connected from below, states ai are connected to above
  UpDownConnectedImaginaryTimePropagation upDownPropagation(betam);
  Wabij.contract(
    1.0, Dai,"ck", Dai,"ai", 0.0,"acik", Bivar_Function<>(upDownPropagation)
  );
  // * interaction V
  Wabij["acik"] *= (*Vabij)["acik"];
  // * thermal weight of contracted indices ck
  Fck["ck"] = 1.0;
  thermalContraction(Fck);
  Wabij["acik"] *= Fck["ck"];
  // T' = W*T
  Tabij["abij"] = Wabij["acik"] * (*amplitudes[0])["cbkj"];
  // finally, freely propagate states bj not connected to V
  FreeImaginaryTimePropagation freePropagation(betam);
  Dai.sum(1.0, Dai,"bj", 0.0,"bj", Univar_Function<>(freePropagation));
  Tabij["abij"] *= Dai["bj"];
  // sum spins * (V*T + T*V) to nextT
  nextT["abij"] =  2.0 * Tabij["abij"];
  nextT["abij"] += 2.0 * Tabij["baji"];

  // T*V*T:
  // effective V = W = propagation of states connected to V
  fetchDelta(Wabij);
  // states cdkl are all connected from below
  SameSideConnectedImaginaryTimePropagation sameSidePropagation(betam);
  Wabij.sum(
    1.0, Wabij,"cdkl", 0.0,"cdkl", Univar_Function<>(sameSidePropagation)
  );
  // * interaction V
  Wabij["cdkl"] *= (*Vabij)["cdkl"];
  // * thermal weight of contracted indices cdkl
  thermalContraction(Wabij);
  Tabij["abij"] =
    (*amplitudes[0])["acik"] * Wabij["cdkl"] * (*amplitudes[0])["dblj"];
  // finally, freely propagate states abij not connected to V
  fetchDelta(Wabij);
  Wabij.sum(1.0, Wabij,"abij", 0.0,"abij", Univar_Function<>(freePropagation));
  Tabij["abij"] *= Wabij["abij"];
  // add spins^2 * T*V*T to nextT
  nextT["abij"] += 4.0 * Tabij["abij"];

  (*amplitudes[recursionLength])["abij"] -= nextT["abij"];
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

