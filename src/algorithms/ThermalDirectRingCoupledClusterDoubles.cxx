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
  Tensor<real> &T0F,
  Tensor<real> &T0FG,
  Tensor<real> &T1F,
  Tensor<real> &T1FG,
  const real DTau,
  Tensor<real> &S1F,
  Tensor<real> &S1FG
) {
  real spins( getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0 );
  ConvolutionC convolutionC(DTau);
  Convolution0 convolution0(DTau);
  Convolution1 convolution1(DTau);
  Convolution00 convolution00(DTau);
  Convolution01 convolution01(DTau);
  Convolution11 convolution11(DTau);

  // TODO: only for real code:
  LOG(1, "FT-DRCCD") << "doubles, constant term..." << std::endl;
  //////////////////////////////////////
  // Sabij = Vabij
  //////////////////////////////////////
  Tensor<real> SFG(*VdFG);
  Transform<real, real>(
    std::function<void(real, real &)>( convolutionC )
  ) (
    (*lambdaFG)["FG"], SFG["FG"]
  );
  S1FG["FG"] -= SFG["FG"];

  LOG(1, "FT-DRCCD") << "doubles, quadratic terms T2 T2..." << std::endl;
  //////////////////////////////////////
  // Sabij = Vklck * Tacik * Tdblj
  //////////////////////////////////////
  //// T^I1(tau_n-1)*T^I2(tau_n-1)
  SFG["FG"] = (+1.0) * spins*spins * T0FG["FH"] * (*VdFG)["HI"] * T0FG["IG"];
  Transform<real, real>(
    std::function<void(real, real &)>( convolution00 )
  ) (
    (*lambdaFG)["FG"], SFG["FG"]
  );
  S1FG["FG"] -= SFG["FG"];
  //// T^I1(tau_n-1)*T^I2(tau_n) and T^I1(tau_n)*T^I2(tau_n-1)
  SFG["FG"] = (+1.0) * spins*spins * T0FG["FH"] * (*VdFG)["HI"] * T1FG["IG"];
  // assemble both T0*T1 and T1*T0
  SFG["FG"] += SFG["GF"];
  Transform<real, real>(
    std::function<void(real, real &)>( convolution01 )
  ) (
    (*lambdaFG)["FG"], SFG["FG"]
  );
  S1FG["FG"] -= SFG["FG"];
  //// T^I1(tau_n)*T^I2(tau_n)
  SFG["FG"] = (+1.0) * spins*spins * T1FG["FH"] * (*VdFG)["HI"] * T1FG["IG"];
  Transform<real, real>(
    std::function<void(real, real &)>( convolution11 )
  ) (
    (*lambdaFG)["FG"], SFG["FG"]
  );
  S1FG["FG"] -= SFG["FG"];
}

