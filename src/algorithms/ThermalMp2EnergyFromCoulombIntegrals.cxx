#include <algorithms/ThermalMp2EnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <math/MultiCombinations.hpp>
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
  beta = 1/getRealArgument("Temperature");

  // compute \Delta^{ab}_{ij} = eps_a + eps_b - eps_i - eps_j
  Dabij = NEW(Tensor<>, false, *Vabij);
  (*Dabij)["abij"] =  (*epsa)["a"];
  (*Dabij)["abij"] += (*epsa)["b"];
  (*Dabij)["abij"] -= (*epsi)["i"];
  (*Dabij)["abij"] -= (*epsi)["j"];

  computeFreeEnergy();
  computeEnergyMoments();
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

void ThermalMp2EnergyFromCoulombIntegrals::computeFreeEnergy() {
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ThermalParticleEigenEnergies"));
  Tensor<> Tabij(*Vabij);
  // apply the 5 terms:
  // Tabij *= f^a = 1/(1+exp(-eps_a*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(ThermalContraction<>(beta, true))
  ) (
    (*epsa)["a"], Tabij["abij"]
  );
  // Tabij *= f^b
  Transform<real, real>(
    std::function<void(real, real &)>(ThermalContraction<>(beta, true))
  ) (
    (*epsa)["b"], Tabij["abij"]
  );
  // Tabij *= f_i = 1/(1+exp(+eps_i*beta))
  Transform<real, real>(
    std::function<void(real, real &)>(ThermalContraction<>(beta, false))
  ) (
    (*epsi)["i"], Tabij["abij"]
  );
  // Tabij *= f_j
  Transform<real, real>(
    std::function<void(real, real &)>(ThermalContraction<>(beta, false))
  ) (
    (*epsi)["j"], Tabij["abij"]
  );
  // Tabij *=
  // integrate(integrate(exp(-Delta*(tau2-tau1),tau2,tau1,beta),tau1,0,beta)
  Transform<real, real>(
    std::function<void(real, real &)>(ThermalMp2Propagation<>(beta))
  ) (
    (*Dabij)["abij"], Tabij["abij"]
  );
  real F(evaluate("F", Tabij, -getRealArgument("Temperature")));
  setRealArgument("ThermalMp2FreeEnergy", F);
}

void ThermalMp2EnergyFromCoulombIntegrals::computeEnergyMoments() {
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  Tensor<> *epsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ThermalParticleEigenEnergies"));

  unsigned int n( std::max(getIntegerArgument("maxEnergyMoment", 2),1l) );
  std::vector<real> energyMoments(n);

  for (unsigned int k(1); k <= n; ++k) {
    // build amplitudes for n-th derivative
    Tensor<> Tabij(false, *Vabij);
    LOG(1, "FT-MP2") << "computing d^" << k << " log(S(beta)) / d(-beta)^" << k
      << " ..." << std::endl;
    // there are 5 terms that can be derived and k derivatives to distribute
    for (auto &derivedTerms: MultiCombinations(5, k)) {
      // convert to degree of derivative for each term
      std::vector<int> degrees(5);
      for (auto term: derivedTerms) { ++degrees[term]; };
      // build contribution from current combination of derived terms
      Tensor<> tabij(*Vabij);
      // apply the derivative of the 5 terms of the respective degrees:
      // Tabij *= f^a = 1/(1+exp(-eps_a*beta))
      Transform<real, real>(
        std::function<void(real, real &)>(
          ThermalContraction<>(beta, true, degrees[0])
        )
      ) (
        (*epsa)["a"], tabij["abij"]
      );
      // Tabij *= f^b
      Transform<real, real>(
        std::function<void(real, real &)>(
          ThermalContraction<>(beta, true, degrees[1])
        )
      ) (
        (*epsa)["b"], tabij["abij"]
      );
      // Tabij *= f_i = 1/(1+exp(+eps_i*beta))
      Transform<real, real>(
        std::function<void(real, real &)>(
          ThermalContraction<>(beta, false, degrees[2])
        )
      ) (
        (*epsi)["i"], tabij["abij"]
      );
      // Tabij *= f_j
      Transform<real, real>(
        std::function<void(real, real &)>(
          ThermalContraction<>(beta, false, degrees[3])
        )
      ) (
        (*epsi)["j"], tabij["abij"]
      );
      // Tabij *=
      // integrate(integrate(exp(-Delta*(tau2-tau1),tau2,tau1,beta),tau1,0,beta)
      Transform<real, real>(
        std::function<void(real, real &)>(
          ThermalMp2Propagation<>(beta, degrees[4])
        )
      ) (
        (*Dabij)["abij"], tabij["abij"]
      );
      // add this contribution
      Tabij["abij"] += tabij["abij"];
    }
    std::stringstream contribution;
    contribution << k << ".central moment of H";
    energyMoments[k-1] = evaluate(contribution.str(), Tabij);
  }
  std::vector<int64_t> indices(Dabij->wrld->rank == 0 ? n : 0);
  for (size_t i(0); i < indices.size(); ++i) { indices[i] = i; };
  auto ctfEnergyMoments(new CTF::Tensor<>(1, std::vector<int>({int(n)}).data()) );
  ctfEnergyMoments->write(indices.size(), indices.data(), energyMoments.data());
  allocatedTensorArgument("ThermalMp2EnergyMoments", ctfEnergyMoments);
}

real ThermalMp2EnergyFromCoulombIntegrals::evaluate(
  const std::string &contribution,
  Tensor<> &Tabij,
  real f
) {
  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Scalar<> energy;
  Tensor<> *Vabij(getTensorArgument("ThermalPPHHCoulombIntegrals"));
  energy[""] = f * (+0.5) * spins * spins * Tabij["abij"] * (*Vabij)["abij"];
  real direct( energy.get_val() );
  energy[""] = f *(-0.5) * spins * Tabij["abij"] * (*Vabij)["abji"];
  real exchange( energy.get_val() );

  LOG(0, "FT-MP2") << contribution << "=" << direct+exchange << std::endl;
  LOG(1, "FT-MP2") << contribution << "_d=" << direct << std::endl;
  LOG(1, "FT-MP2") << contribution << "_x=" << exchange << std::endl;

  return direct+exchange;
}

