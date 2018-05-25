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
  CTF::Tensor<real> &T0FG,
  CTF::Tensor<real> &T1FG,
  const real DTau,
  CTF::Tensor<real> &S1FG
) {
  // TODO: only for real code:
  // constant term:
  real spins(2.0);
  LOG(1, "FT-DRCCD") << "constant term..." << std::endl;
  CTF::Tensor<real> SFG(*VdFG);
  Transform<real, real>(
    std::function<void(real, real &)>( ConvolutionC(DTau) )
  ) (
    (*lambdaFG)["FG"], SFG["FG"]
  );
  S1FG["FG"] -= SFG["FG"];

  LOG(1, "FT-DRCCD") << "quadratic terms..." << std::endl;
  // quadratic terms:
  //   T^I1(tau_n-1)*T^I2(tau_n-1)
  SFG["FG"] = T0FG["FH"] * (*VdFG)["HI"] * T0FG["IG"];
  SFG["FG"] *= spins*spins;
  Transform<real, real>(
    std::function<void(real, real &)>( Convolution00(DTau) )
  ) (
    (*lambdaFG)["FG"], SFG["FG"]
  );
  S1FG["FG"] -= SFG["FG"];
  //   T^I1(tau_n-1)*T^I2(tau_n) and T^I1(tau_n)*T^I2(tau_n-1)
  SFG["FG"] =  T0FG["FH"] * (*VdFG)["HI"] * T1FG["IG"];
  SFG["FG"] += T1FG["FH"] * (*VdFG)["HI"] * T0FG["IG"];
  SFG["FG"] *= spins*spins;
  Transform<real, real>(
    std::function<void(real, real &)>( Convolution01(DTau) )
  ) (
    (*lambdaFG)["FG"], SFG["FG"]
  );
  S1FG["FG"] -= SFG["FG"];
  //   T^I1(tau_n)*T^I2(tau_n)
  SFG["FG"] =  T1FG["FH"] * (*VdFG)["HI"] * T1FG["IG"];
  SFG["FG"] *= spins*spins;
  Transform<real, real>(
    std::function<void(real, real &)>( Convolution11(DTau) )
  ) (
    (*lambdaFG)["FG"], SFG["FG"]
  );
  S1FG["FG"] -= SFG["FG"];

  // shift back
//  S1FG["FG"] += 0.1;
}

