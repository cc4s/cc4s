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
  Tensor<real> &T0abij,
  Tensor<real> &T1abij,
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

  // TODO: only for real code:
  LOG(1, "FT-DRCCD") << "doubles, constant term..." << std::endl;
  //////////////////////////////////////
  // Sabij = Vabij
  //////////////////////////////////////
  Tensor<real> SFG(*VdFG);
  Transform<real, real>(std::function<void(real, real &)>(convolutionC))(
    (*lambdaFG)["FG"], SFG["FG"]
  );
  S1FG["FG"] -= SFG["FG"];

  LOG(1, "FT-DRCCD") << "doubles, quadratic terms T2 T2..." << std::endl;
  //////////////////////////////////////
  // Sabij = Vklck * Tacik * Tdblj
  //////////////////////////////////////
  auto Vabij( getTensorArgument<real>("ThermalPPHHCoulombIntegrals") );
  Tensor<real> Sabij(false, *Vabij);
  //// T^I1(tau_n-1)*T^I2(tau_n-1)
  Sabij["abij"] = (+1.0) * spins*spins *
    T0abij["acik"] * (*Vabij)["cdkl"] * T0abij["dblj"];
  propagateAmplitudes(Sabij, convolution00, S1FG);

  //// T^I1(tau_n-1)*T^I2(tau_n) and T^I1(tau_n)*T^I2(tau_n-1)
  Sabij["abij"] = (+1.0) * spins*spins *
    T0abij["acik"] * (*Vabij)["cdkl"] * T1abij["dblj"];
  // assemble both T0*T1 and T1*T0
  Sabij["abij"] += Sabij["baji"];
  propagateAmplitudes(Sabij, convolution01, S1FG);

  //// T^I1(tau_n)*T^I2(tau_n)
  Sabij["abij"] = (+1.0) * spins*spins *
    T1abij["acik"] * (*Vabij)["cdkl"] * T1abij["dblj"];
  propagateAmplitudes(Sabij, convolution11, S1FG);
}

