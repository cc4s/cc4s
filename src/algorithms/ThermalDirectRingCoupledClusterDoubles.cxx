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
  Tensor<real> &T0FG,
  Tensor<real> &T1FG,
  const real DTau,
  Tensor<real> &S1FG
) {
  real spins( getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0 );
  ConvolutionC convolutionC(DTau);
  Convolution0 convolution0(DTau);
  Convolution1 convolution1(DTau);
  Convolution00 convolution00(DTau);
  Convolution01 convolution01(DTau);
  Convolution11 convolution11(DTau);
/*
  auto Vabij( getTensorArgument<real>("ThermalPPHHCoulombIntegrals") );
  auto Na( getTensorArgument<real>("ThermalParticleOccupancies") );
  auto Ni( getTensorArgument<real>("ThermalHoleOccupancies") );
*/
  Transform<real> chop(
    std::function<void(real &)>(
      [](real &t) {
        if (std::abs(t) < 64*sqrt(std::numeric_limits<real>::epsilon())) t = 0.0;
      }
    )
  );

  // TODO: only for real code:
  LOG(1, "FT-DRCCD") << "doubles, constant term..." << std::endl;
  //////////////////////////////////////
  // Sabij = Vabij
  //////////////////////////////////////
  Tensor<real> SFG(*VdFG);
//    (*ga)["a"] * (*ga)["b"] * (*gi)["i"] * (*gi)["j"];
//    (*Na)["a"] * (*Na)["b"] * (*Ni)["i"] * (*Ni)["j"];
//  chop(Sabij["abij"]);
  propagateAmplitudes(SFG, convolutionC);
  S1FG["FG"] += (-1.0) * SFG["FG"];

//  return;
/*
  if (isArgumentGiven("ThermalPHPHCoulombIntegrals")) {
    // particle/hole ladder of H1, if given
    auto Vbiaj(getTensorArgument<>("ThermalPHPHCoulombIntegrals"));
    Sabij["abij"] = (-1.0) * T0abij["ackj"] * (*Vbiaj)["bkci"];
    chop(Sabij["abij"]);
    Sabij["abij"] += Sabij["baji"];
    propagateAmplitudes(Sabij, convolution0);
    S1abij["abij"] += (-1.0) * Sabij["abij"];

    Sabij["abij"] = (-1.0) * T1abij["ackj"] * (*Vbiaj)["bkci"];
    chop(Sabij["abij"]);
    Sabij["abij"] += Sabij["baji"];
    propagateAmplitudes(Sabij, convolution1);
    S1abij["abij"] += (-1.0) * Sabij["abij"];
  }
*/

  LOG(1, "FT-DRCCD") << "doubles, quadratic terms T2 T2..." << std::endl;
  //////////////////////////////////////
  // Sabij = Vklck * Tacik * Tdblj
  //////////////////////////////////////
  Tensor<real> WFG(false, *VdFG);
  WFG["FG"]  = (+1.0) * spins*spins * (*VdFG)["FG"];
  if (getIntegerArgument("adjacentPairsExchange", 0)) {
    WFG["FG"] += (-1.0) * spins * (*VxFG)["FG"];
  }
  //// T^I1(tau_n-1)*T^I2(tau_n-1)
  SFG["FG"] = T0FG["FH"] * WFG["HI"] * T0FG["IG"];
  propagateAmplitudes(SFG, convolution00);
  chop(SFG["FG"]);
  S1FG["FG"] += (-1.0) * SFG["FG"];

  //// T^I1(tau_n-1)*T^I2(tau_n) and T^I1(tau_n)*T^I2(tau_n-1)
  SFG["FG"] = T0FG["FH"] * WFG["HI"] * T1FG["IG"];
  // assemble both T0*T1 and T1*T0
  chop(SFG["FG"]);
  SFG["FG"] += SFG["GF"];
  propagateAmplitudes(SFG, convolution01);
  S1FG["FG"] += (-1.0) * SFG["FG"];

  //// T^I1(tau_n)*T^I2(tau_n)
  SFG["FG"] = T1FG["FH"] * WFG["HI"] * T1FG["IG"];
  chop(SFG["FG"]);
  propagateAmplitudes(SFG, convolution11);
  S1FG["FG"] += (-1.0) * SFG["FG"];
}

