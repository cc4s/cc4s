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
  int vo[] = {Nv, No};
  auto fVaijb(NEW(Tensor<real>, 4, voov, Vabij->sym));
  auto Vijab(NEW(Tensor<real>, 4, oovv, Vabij->sym));
  CTF::Tensor<real> fai(2, vo, &Vabij->sym[1]);
  thermalContraction(fai);
  (*fVaijb)["ajib"] = (*Vabij)["abij"] * fai["ai"];
  (*Vijab)["ijab"] = (*Vabij)["abij"];

  // constant term:
  real spins(2.0);
  LOG(1, "FT-DRCCD") << "constant term..." << std::endl;  
  CTF::Tensor<real> Sabij(*Vabij);
  thermalContraction(Sabij);
  Transform<real, real>(
    std::function<void(real, real &)>( ConvolutionC(DTau) )
  ) (
    (*Dabij)["abij"], Sabij["abij"]
  );
  S1abij["abij"] -= Sabij["abij"];

/*
  // linear terms:
  LOG(1, "FT-DRCCD") << "linear terms..." << std::endl;  
  //   T^I(tau_n-1):
  Sabij["abij"] =  T0abij["acik"] * (*fVaijb)["bkjc"];
  Sabij["abij"] += T0abij["dblj"] * (*fVaijb)["alid"];
  Sabij["abij"] *= spins;
  Transform<real, real>(
    std::function<void(real, real &)>( Convolution0(DTau) )
  ) (
    (*Dabij)["abij"], Sabij["abij"]
  );
  S1abij["abij"] -= Sabij["abij"];
  //   T^I(tau_n):
  Sabij["abij"] =  T1abij["acik"] * (*fVaijb)["bkjc"];
  Sabij["abij"] += T1abij["dblj"] * (*fVaijb)["alid"];
  Sabij["abij"] *= spins;
  Transform<real, real>(
    std::function<void(real, real &)>( Convolution1(DTau) )
  ) (
    (*Dabij)["abij"], Sabij["abij"]
  );
  S1abij["abij"] -= Sabij["abij"];

  LOG(1, "FT-DRCCD") << "quadratic terms..." << std::endl;  

  // quadratic terms:
  //   T^I1(tau_n-1)*T^I2(tau_n-1)
  Sabij["abij"] = T0abij["acik"] * (*Vijab)["klcd"] * T0abij["dblj"];
  Sabij["abij"] *= spins*spins;
  Transform<real, real>(
    std::function<void(real, real &)>( Convolution00(DTau) )
  ) (
    (*Dabij)["abij"], Sabij["abij"]
  );
  S1abij["abij"] -= Sabij["abij"];
  //   T^I1(tau_n-1)*T^I2(tau_n) and T^I1(tau_n)*T^I2(tau_n-1)
  Sabij["abij"] =  T0abij["acik"] * (*Vijab)["klcd"] * T1abij["dblj"];
  Sabij["abij"] += T1abij["acik"] * (*Vijab)["klcd"] * T0abij["dblj"];
  Sabij["abij"] *= spins*spins;
  Transform<real, real>(
    std::function<void(real, real &)>( Convolution01(DTau) )
  ) (
    (*Dabij)["abij"], Sabij["abij"]
  );
  S1abij["abij"] -= Sabij["abij"];
  //   T^I1(tau_n)*T^I2(tau_n)
  Sabij["abij"] =  T1abij["acik"] * (*Vijab)["klcd"] * T1abij["dblj"];
  Sabij["abij"] *= spins*spins;
  Transform<real, real>(
    std::function<void(real, real &)>( Convolution11(DTau) )
  ) (
    (*Dabij)["abij"], Sabij["abij"]
  );
*/
}

