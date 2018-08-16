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
  auto epsi( getTensorArgument<>("ThermalHoleEigenEnergies") );
  auto epsa( getTensorArgument<>("ThermalParticleEigenEnergies") );
  auto deltaai( getTensorArgument<>("ThermalParticleHoleOverlap") );
  int Nv(epsa->lens[0]), No(epsi->lens[0]), NF(S1F.lens[0]);
  Tensor<real> T0ai(2, std::vector<int>({Nv,No}).data());
  Tensor<real> T1ai(T0ai);
  Tensor<real> H0ai(T0ai);
  T0ai["ai"] = (*UaiF)["aiF"] * T0F["F"];
  T1ai["ai"] = (*UaiF)["aiF"] * T1F["F"];
  H0ai["ai"] = (*gi)["i"] * (*deltaai)["ai"] * (*ga)["a"] * (*epsa)["a"];

  // TODO: only for real code:
  //////////////////////////////////////////////////////////////////////////////
  // singles:
  //////////////////////////////////////////////////////////////////////////////
  LOG(1, "FT-DRCCD") << "singles, constant term..." << std::endl;
  //////////////////////////////////////
  // Sai = H0ai
  //////////////////////////////////////
  Tensor<real> SF(*H0F);
  Transform<real, real>(
    std::function<void(real, real &)>( convolutionC )
  ) (
    (*lambdaF)["F"], SF["F"]
  );
  S1F["F"] -= SF["F"];

/*
  LOG(1, "FT-DRCCD") << "singles, linear terms..." << std::endl;
  SF["F"] = (+1.0) * spins * (*H0F)["G"] * T0FG["FG"];
  Transform<real, real>(
    std::function<void(real, real &)>( convolution0 )
  ) (
    (*lambdaF)["F"], SF["F"]
  );
  S1F["F"] -= SF["F"];
  SF["F"] = 2 * (*H0F)["G"] * T1FG["FG"];
  Transform<real, real>(
    std::function<void(real, real &)>( convolution1 )
  ) (
    (*lambdaF)["F"], SF["F"]
  );
  S1F["F"] -= SF["F"];
*/
  LOG(1, "FT-DRCCD") << "singles, quadratic terms..." << std::endl;
  //////////////////////////////////////
  // Sai = H0kc Tak Tci
  //////////////////////////////////////
  //// T^I1(tau_n-1)*T^I2(tau_n-1)
  SF["F"] = (-1.0) * (*UaiF)["aiF"] * T0ai["ak"] * H0ai["ck"] * T0ai["ci"];
  Transform<real, real>(
    std::function<void(real, real &)>( convolution00 )
  ) (
    (*lambdaF)["F"], SF["F"]
  );
  S1F["F"] -= SF["F"];
  //// T^I1(tau_n-1)*T^I2(tau_n) and T^I1(tau_n)*T^I2(tau_n-1)
  SF["F"]  = (-1.0) * (*UaiF)["aiF"] * T0ai["ak"] * H0ai["ck"] * T1ai["ci"];
  SF["F"] += (-1.0) * (*UaiF)["aiF"] * T1ai["ak"] * H0ai["ck"] * T0ai["ci"];
  Transform<real, real>(
    std::function<void(real, real &)>( convolution01 )
  ) (
    (*lambdaF)["F"], SF["F"]
  );
  S1F["F"] -= SF["F"];
  //// T^I1(tau_n)*T^I2(tau_n)
  SF["F"] = (-1.0) * (*UaiF)["aiF"] * T1ai["ak"] * H0ai["ck"] * T1ai["ci"];
  Transform<real, real>(
    std::function<void(real, real &)>( convolution11 )
  ) (
    (*lambdaF)["F"], SF["F"]
  );
  S1F["F"] -= SF["F"];

  //////////////////////////////////////////////////////////////////////////////
  // doubles:
  //////////////////////////////////////////////////////////////////////////////
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

  LOG(1, "FT-DRCCD") << "doubles, quadratic terms T1 T2..." << std::endl;
  //////////////////////////////////////
  // Sabij = Hkc * Tak * Tcbij
  //////////////////////////////////////
  Tensor<real> T0aiG(3, std::vector<int>({Nv,No,NF}).data());
  Tensor<real> T1aiG(T0aiG);
  T0aiG["aiG"] = (*UaiF)["aiF"] * T0FG["FG"];
  T1aiG["aiG"] = (*UaiF)["aiF"] * T1FG["FG"];
  //// T^I1(tau_n-1)*T^I2(tau_n-1)
  SFG["FG"] = (-1.0) * (*UaiF)["aiF"] * T0ai["ak"] * H0ai["ck"] * T0aiG["ciG"];
  // symmetrize left/right
  SFG["FG"] += SFG["GF"];
  Transform<real, real>(
    std::function<void(real, real &)>( convolution00 )
  ) (
    (*lambdaFG)["FG"], SFG["FG"]
  );
  S1FG["FG"] -= SFG["FG"];
  //// T^I1(tau_n-1)*T^I2(tau_n) and T^I1(tau_n)*T^I2(tau_n-1)
  SFG["FG"]  = (-1.0) * (*UaiF)["aiF"] * T0ai["ak"] * H0ai["ck"] * T1aiG["ciG"];
  SFG["FG"] += (-1.0) * (*UaiF)["aiF"] * T1ai["ak"] * H0ai["ck"] * T0aiG["ciG"];
  // symmetrize left/right
  SFG["FG"] += SFG["GF"];
  Transform<real, real>(
    std::function<void(real, real &)>( convolution01 )
  ) (
    (*lambdaFG)["FG"], SFG["FG"]
  );
  S1FG["FG"] -= SFG["FG"];
  //// T^I1(tau_n)*T^I2(tau_n)
  SFG["FG"] = (-1.0) * (*UaiF)["aiF"] * T1ai["ak"] * H0ai["ck"] * T1aiG["ciG"];
  // symmetrize left/right
  SFG["FG"] += SFG["GF"];
  Transform<real, real>(
    std::function<void(real, real &)>( convolution11 )
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

