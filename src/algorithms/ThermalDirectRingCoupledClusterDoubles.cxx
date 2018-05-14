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
  CTF::Tensor<real> &T0abij,
  CTF::Tensor<real> &T1abij,
  const real DTau,
  CTF::Tensor<real> &S1abij
) {
  auto Vabij(getTensorArgument<real>("ThermalPPHHCoulombIntegrals"));
  // TODO: only for real code:
  int Nv(Vabij->lens[0]), No(Vabij->lens[2]);
  int voov[] = {Nv, No, No, Nv};
  int oovv[] = {No, No, Nv, Nv};
  auto Vaijb(NEW(Tensor<real>, 4, voov, Vabij->sym));
  auto Vijab(NEW(Tensor<real>, 4, oovv, Vabij->sym));
  (*Vaijb)["ajib"] = (*Vabij)["abij"];
  (*Vijab)["ijab"] = (*Vabij)["abij"];

  // constant term:
  CTF::Tensor<real> Sabij(*Vabij);
  thermalContraction(Sabij);
  Transform<real, real>(
    std::function<void(real, real &)>( ConvolutionC(DTau) )
  ) (
    (*Dabij)["abij"], Sabij["abij"]
  );
  S1abij["abij"] -= Sabij["abij"];

  // linear terms:
  CTF::Tensor<real> fai(2, &Vabij->lens[1], &Vabij->sym[1]);
  thermalContraction(fai);
  //   T^I(tau_n-1):
  Sabij["abij"] =  T0abij["acik"] * (*Vaijb)["bkjc"] * fai["bj"];
  Sabij["abij"] += T0abij["dblj"] * (*Vaijb)["alid"] * fai["ai"];
  Transform<real, real>(
    std::function<void(real, real &)>( Convolution0(DTau) )
  ) (
    (*Dabij)["abij"], Sabij["abij"]
  );
  S1abij["abij"] -= Sabij["abij"];
  //   T^I(tau_n):
  Sabij["abij"] =  T1abij["acik"] * (*Vaijb)["bkjc"] * fai["bj"];
  Sabij["abij"] += T1abij["dblj"] * (*Vaijb)["alid"] * fai["ai"];
  Transform<real, real>(
    std::function<void(real, real &)>( Convolution1(DTau) )
  ) (
    (*Dabij)["abij"], Sabij["abij"]
  );
  S1abij["abij"] -= Sabij["abij"];
  
  // quadratic terms:
  //   T^I1(tau_n-1)*T^I2(tau_n-1)
  Sabij["abij"] =  T0abij["acik"] * (*Vijab)["klcd"] * T0abij["dblj"];
  Transform<real, real>(
    std::function<void(real, real &)>( Convolution00(DTau) )
  ) (
    (*Dabij)["abij"], Sabij["abij"]
  );
  S1abij["abij"] -= Sabij["abij"];
  //   T^I1(tau_n-1)*T^I2(tau_n) and T^I1(tau_n)*T^I2(tau_n-1)
  Sabij["abij"] =  T0abij["acik"] * (*Vijab)["klcd"] * T1abij["dblj"];
  Sabij["abij"] += T1abij["acik"] * (*Vijab)["klcd"] * T0abij["dblj"];
  Transform<real, real>(
    std::function<void(real, real &)>( Convolution01(DTau) )
  ) (
    (*Dabij)["abij"], Sabij["abij"]
  );
  S1abij["abij"] -= Sabij["abij"];
  //   T^I1(tau_n)*T^I2(tau_n)
  Sabij["abij"] =  T1abij["acik"] * (*Vijab)["klcd"] * T1abij["dblj"];
  Transform<real, real>(
    std::function<void(real, real &)>( Convolution11(DTau) )
  ) (
    (*Dabij)["abij"], Sabij["abij"]
  );
}

