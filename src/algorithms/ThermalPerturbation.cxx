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
  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ThermalParticleEigenEnergies"));
  Tensor<> *Vaibj(getTensorArgument("ThermalPHPHCoulombIntegrals"));
  Tensor<> *Vaijk(getTensorArgument("ThermalPHHHCoulombIntegrals"));
  Tensor<> *Vijkl(getTensorArgument("ThermalHHHHCoulombIntegrals"));
  int No(epsi->lens[0]), Nv(epsa->lens[0]);
  beta = 1 / getRealArgument("Temperature");
  deltaMu = getRealArgument("chemicalPotentialShift", 0.0);

  // compute contraction weights
  Tensor<> Nc(1, &Nv);
  Tensor<> Nk(1, &No);
  Tensor<> nk(1, &No);
  // particles in perturbation: contraction with shifted chemical potential
  Nc["c"] = 1.0;
  // Nc *= f^c = 1/(1+exp(-(eps_c-deltaMu)*beta))
  // NOTE: eps_q already contains epsilon_q - mu_0,
  // where mu_0 is the Hartree--Fock chemical potential
  Transform<real, real>(
    std::function<void(real, real &)>(
      ThermalContraction<>(beta, deltaMu, true)
    )
  ) (
    (*epsa)["c"], Nc["c"]
  );
  (*getTensorArgument("ThermalParticleOccupancies"))["c"] = Nc["c"];
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
  (*getTensorArgument("ThermalHoleOccupancies"))["k"] = Nk["k"];
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
  real Omega0(-getDLogZH0(0,0)/beta);

  // first order:
  // Hartree and exchange term, use shifted occupancies Nk
  // both terms have left/right symmetry
  Scalar<> energy;
  energy[""] = (+0.5) * spins * spins * Nk["i"] * Nk["j"] * (*Vijkl)["ijij"];
  real Omega1D( energy.get_val() );
  energy[""] = (-0.5) * spins * Nk["i"] * Nk["j"] * (*Vijkl)["ijji"];
  real Omega1X( energy.get_val() );
  // minus effective potential,
  // use Hartree--Fock occupancies nk for the contraction k, composing Veff
  // use occupancies Nk with shifted mu for contraction i, using Veff
  energy[""] = (-1.0) * spins * spins * Nk["i"] * nk["k"] * (*Vijkl)["ikik"];
  real Omega1Deff( energy.get_val() );
  energy[""] = (+1.0) * spins * Nk["i"] * nk["k"] * (*Vijkl)["ikki"];
  real Omega1Xeff( energy.get_val() );

  // one-body part of perturbation for higher orders:
  // contraction weight = weight of perturation - weight of effective pot.
  Nk["k"] += (-1.0) * nk["k"];
  // hole-hole:
  Tensor<> *Vij(new Tensor<>(2, std::vector<int>({No,No}).data()));
  allocatedTensorArgument<>("ThermalHHPerturbation", Vij);
  (*Vij)["ij"] =  (+1.0) * spins * (*Vijkl)["ikjk"] * Nk["k"];
  (*Vij)["ij"] += (-1.0) * (*Vijkl)["ikkj"] * Nk["k"];
  // particle-hole:
  Tensor<> *Vai(new Tensor<>(2, std::vector<int>({Nv,No}).data()));
  allocatedTensorArgument<>("ThermalPHPerturbation", Vai);
  (*Vai)["ai"] =  (+1.0) * spins * (*Vaijk)["akik"] * Nk["k"];
  (*Vai)["ai"] += (-1.0) * (*Vaijk)["akki"] * Nk["k"];
  // particle-particle:
  Tensor<> *Vab(new Tensor<>(2, std::vector<int>({Nv,Nv}).data()));
  allocatedTensorArgument<>("ThermalPPPerturbation", Vab);
  (*Vab)["ab"] =  (+1.0) * spins * (*Vaibj)["akbk"] * Nk["k"];
  (*Vab)["ab"] += (-1.0) * (*Vaijk)["akki"] * Nk["k"];

  real OmegaHf(Omega0+Omega1D+Omega1X+Omega1Deff+Omega1Xeff);
  EMIT() << YAML::Key << "Omega0" << YAML::Value << Omega0;
  EMIT() << YAML::Key << "Omega1D" << YAML::Value << Omega1D;
  EMIT() << YAML::Key << "Omega1X" << YAML::Value << Omega1X;
  EMIT() << YAML::Key << "Omega1Deff" << YAML::Value << Omega1Deff;
  EMIT() << YAML::Key << "Omega1Xeff" << YAML::Value << Omega1Xeff;
  EMIT() << YAML::Key << "Hartree-Fock-grand-potential"
    << YAML::Value << OmegaHf;
}

void ThermalPerturbation::dryRun() {
  // FIXME
}


cc4s::real ThermalPerturbation::getDLogZH0(
  const unsigned int dbeta_n, const unsigned int dmu_m
) {
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> Ti(false, *epsi);
  Ti["i"] = 1;
  switch (dbeta_n+dmu_m) {
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
    // FIXME: only works for dbeta_n+dmu_m=1
    Transform<real, real>(
      std::function<void(real, real &)>(
        ThermalContraction<>(beta, deltaMu, false, 0, 0)
      )
    ) (
      (*epsi)["i"], Ti["i"]
    );
    if (dbeta_n>0) {
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

