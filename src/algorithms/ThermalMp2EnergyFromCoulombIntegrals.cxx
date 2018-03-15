#include <algorithms/ThermalMp2EnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ThermalMp2EnergyFromCoulombIntegrals);

ThermalMp2EnergyFromCoulombIntegrals::ThermalMp2EnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ThermalMp2EnergyFromCoulombIntegrals::~ThermalMp2EnergyFromCoulombIntegrals() {
}

void ThermalMp2EnergyFromCoulombIntegrals::run() {
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ThermalParticleEigenEnergies"));
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  real beta( 1/getRealArgument("Temperature") );

  // first, compute \Delta^{ab}_{ij}, note that its sign is opposite than usual
  Tensor<> Dabij(false, *Vabij);
  Dabij["abij"] =  (*epsa)["a"];
  Dabij["abij"] += (*epsa)["b"];
  Dabij["abij"] -= (*epsi)["i"];
  Dabij["abij"] -= (*epsi)["j"];

  Tensor<> Tabij(*Vabij);
  Transform<real, real>(
    std::function<void(real, real &)>(ThermalContraction<>(+beta))
  ) (
    (*epsa)["a"], Tabij["abij"]
  );
  Transform<real, real>(
    std::function<void(real, real &)>(ThermalContraction<>(+beta))
  ) (
    (*epsa)["b"], Tabij["abij"]
  );
  Transform<real, real>(
    std::function<void(real, real &)>(ThermalContraction<>(-beta))
  ) (
    (*epsi)["i"], Tabij["abij"]
  );
  Transform<real, real>(
    std::function<void(real, real &)>(ThermalContraction<>(-beta))
  ) (
    (*epsi)["j"], Tabij["abij"]
  );
  Transform<real, real>(
    std::function<void(real, real &)>(ThermalMp2Propagation<>(beta))
  ) (
    Dabij["abij"], Tabij["abij"]
  );

  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Scalar<> energy;
  energy[""] = +0.5 * spins * spins * Tabij["abij"] * (*Vabij)["abij"];
  real direct( -energy.get_val()/beta );
  energy[""] = -0.5 * spins * Tabij["abij"] * (*Vabij)["abji"];
  real exchange( -energy.get_val()/beta );

  LOG(0, "FT-MP2") << "FT-MP2=" << direct+exchange << std::endl;
  LOG(1, "FT-MP2") << "FT-MP2d=" << direct << std::endl;
  LOG(1, "FT-MP2") << "FT-MP2x=" << exchange << std::endl;

  setRealArgument("ThermalMp2Energy", direct+exchange);
}

void ThermalMp2EnergyFromCoulombIntegrals::dryRun() {
  //DryTensor<> *Vabij(
  getTensorArgument<double, DryTensor<double>>("ThermalPPHHCoulombIntegrals");
  //);

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("ThermalHoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ThermalParticleEigenEnergies")
  );
  
  // Compute the No,Nv
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Allocate the doubles amplitudes
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { Nv, Nv, No, No };
  DryTensor<> Tabij(4, vvoo, syms);

  DryScalar<> energy();
}

