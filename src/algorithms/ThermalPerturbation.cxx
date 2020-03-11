/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <algorithms/ThermalPerturbation.hpp>
#include <math/MathFunctions.hpp>
#include <math/MultiCombinations.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Emitter.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ThermalPerturbation);

ThermalPerturbation::ThermalPerturbation(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ThermalPerturbation::~ThermalPerturbation() {
}

void ThermalPerturbation::run() {
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ThermalParticleEigenEnergies"));
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  beta = 1 / getRealArgument("Temperature");

  // compute \Delta^{ab}_{ij} = eps_a + eps_b - eps_i - eps_j
  Dabij = NEW(Tensor<>, false, *Vabij);
  (*Dabij)["abij"] =  (*epsa)["a"];
  (*Dabij)["abij"] += (*epsa)["b"];
  (*Dabij)["abij"] -= (*epsi)["i"];
  (*Dabij)["abij"] -= (*epsi)["j"];
  Dai = NEW(Tensor<>, 2, &Vabij->lens[1]);
  (*Dai)["ai"] =  (*epsa)["a"];
  (*Dai)["ai"] -= (*epsi)["i"];

  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Tensor<> Nc(1, epsa->lens);
  Tensor<> Nk(1, epsi->lens);
  Tensor<> nk(1, epsi->lens);
  // particles in perturbation: contraction with shifted chemical potential
  Nc["c"] = 1.0;
  // Nc *= f^c = 1/(1+exp(-(eps_c-deltaMu)*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, deltaMu, true)
    )
  ) (
    (*epsa)["c"], Nc["c"]
  );
  // holes in perturbation: contraction with shifted chemical potential
  Nk["k"] = 1.0;
  // Nk *= f_k = 1/(1+exp(+(eps_k-deltaMu)*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, deltaMu, false)
    )
  ) (
    (*epsi)["k"], Nk["k"]
  );
  // terms in effective potential: contraction with Hartree--Fock chemical pot.
  nk["k"] = 1.0;
  // nk *= f_k = 1/(1+exp(+(eps_k-0)*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, 0.0, false)
    )
  ) (
    (*epsi)["k"], nk["k"]
  );

  // zeroth order
  real Omega0(-getDLogZH0(0, D_BETA)/beta);

  // first order
  // Hartree and exchange term, use shifted occupancies Nk
  Tensor<> *Vijkl(getTensorArgument("ThermalHHHHCoulombIntegrals"));
  energy[""] = (+0.5) * spins * spins * Nk["i"] * Nk["j"] * (*Vijkl)["ijij"];
  real ED1( energy.get_val() );
  energy[""] = (-0.5) * spins * Nk["i"] * Nk["j"] * (*Vijkl)["ijji"];
  real EX1( energy.get_val() );
  // minus effective potential, use Hartree--Fock occupancies nk
  energy[""] = (-1.0) * spins * spins * nk["i"] * nk["j"] * (*Vijkl)["ijij"];
  real EE1( energy.get_val() );
  energy[""] = (+1.0) * spins * Nk["i"] * Nk["j"] * (*Vijkl)["ijji"];

  // second order:
  Tensor<> *Vaijk(getTensorArgument("ThermalPHHHCoulombIntegrals"));
  // operator for singles
  Tensor<> Fai(2, Vaijk->lens);
  // contraction weight = weight of perturation - weight of effective pot.
  Nk["k"] -= nk["k"];
  // setup singles operator with difference of perturbation - effective pot.
  Fai["ai"] =  (+1.0) * spins * (*Vaijk)["akik"] * Nk["k"];
  Fai["ai"] += (-1.0) * (*Vaijk)["akki"] * Nk["k"];
  // restore normal weight of pertrubation holes for further calculations
  Nk["k"] += nk["k"];
  // singles:
  // start with Fai
  Tensor<> Tai(Fai);
  // Tai *=
  // integrate(integrate(exp(-Delta*(tau2-tau1),tau2,tau1,beta),tau1,0,beta)
  Transform<real, real>(
    std::function<void(real, real &)>(ThermalMp2Propagation<>(beta, 0))
  ) (
    (*Dai)["ai"], Tai["ai"]
  );
  // no symmetry, one loop, one hole contracted: +1.0 * spins
  energy[""] = (+1.0) * spins * Tai["ai"] * Fai["ai"] * Nk["i"] * Nc["c"];
  real ES2( -energy.get_val()/beta );


  real FHf(Omega0+ED1+EX1+EE1);
  EMIT() << YAML::Key << "Omega0" << YAML::Value << Omega0;
  EMIT() << YAML::Key << "D1" << YAML::Value << ED1;
  EMIT() << YAML::Key << "X1" << YAML::Value << EX1;
  EMIT() << YAML::Key << "eff1" << YAML::Value << EE1;
  EMIT() << YAML::Key << "Hartree-Fock-grand-potential" << YAML::Value << FHf;
}

void ThermalPerturbation::dryRun() {
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
}


cc4s::real ThermalPerturbation::getDLogZH0(
  const unsigned int n, const bool dbeta
) {
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> Ti(false, *epsi);
  Ti["i"] = 1;
  switch (n) {
  case 0:
    Transform<real, real>(
      std::function<void(real, real &)>(
        [this](real eps, real &T) {
          T = std::log(1 + std::exp(-beta*(eps-deltaMu)));
        }
      )
    ) (
      (*epsi)["i"], Ti["i"]
    );
    break;
  default:
    // the first derivative is eps_i*f_i (dbeta) or f_i (dmu),
    // use the ThermalContraction clas providing arbitrary derivatives of f_i
    Transform<real, real>(
      std::function<void(real, real &)>(
        ThermalContraction<>(beta, deltaMu, false, n-1, dbeta)
      )
    ) (
      (*epsi)["i"], Ti["i"]
    );
    if (dbeta) {
      Ti["i"] *= (*epsi)["i"];
    }
    break;
  }
  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Scalar<> dLogZH0;
  // sum over all relevant (thermal hole) states
  dLogZH0[""] = spins * Ti["i"];
  return dLogZH0.get_val();
}

