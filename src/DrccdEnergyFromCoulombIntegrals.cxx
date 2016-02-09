#include <DrccdEnergyFromCoulombIntegrals.hpp>
#include <util/Log.hpp>
#include <util/MathFunctions.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <Options.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DrccdEnergyFromCoulombIntegrals);

DrccdEnergyFromCoulombIntegrals::DrccdEnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

DrccdEnergyFromCoulombIntegrals::~DrccdEnergyFromCoulombIntegrals() {
}

/**
 * \brief Calculates MP2 energy from Coulomb integrals Vabij
 */
void DrccdEnergyFromCoulombIntegrals::run() {
  Tensor<> *Vabij(getTensorArgument("ParticleHoleCoulombIntegrals"));

  int nv(Vabij->lens[0]);
  int no(Vabij->lens[2]);
  int lens[] = { nv, nv, no, no };
  int syms[] = { NS, NS, NS, NS };
  Tensor<> *Tabij(new Tensor<>(4, lens, syms, *Cc4s::world, "Tabij"));
  allocatedTensorArgument("DrccdDoublesAmplitudes", Tabij);
  
  Scalar<> energy(*Cc4s::world);
  double e(0), dire, exce;

  LOG(0) <<
    "Solving direct ring Coupled Cluster Doubles Amplitude Equations:" <<
    std::endl;
 
  for (int i(0); i < Cc4s::options->niter; ++i) {
    LOG(0) << "iteration: " << i << std::endl;
    iterate();
    energy[""] = 2.0 * (*Tabij)["abij"] * (*Vabij)["abij"];
    dire = energy.get_val();
    energy[""] = (*Tabij)["abji"] * (*Vabij)["abij"];
    exce = -1.0 * energy.get_val();
    e = dire + exce;
    LOG(0) << "e=" << e << std::endl;
    LOG(1) << "RPA=" << dire << std::endl;
    LOG(1) << "SOSEX=" << exce << std::endl;
  }

  setRealArgument("DrccdEnergy", e);
}

void DrccdEnergyFromCoulombIntegrals::iterate() {
  // get tensors
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<> *Vabij(getTensorArgument("ParticleHoleCoulombIntegrals"));
  Tensor<> *Tabij(getTensorArgument("DrccdDoublesAmplitudes"));
  Tensor<> Rabij(false, *Vabij);
  Tensor<> Cabij(false, *Vabij);
  Tensor<> Dabij(false, *Vabij);

  Rabij["abij"] = (*Vabij)["abij"];
  Rabij["abij"] += 2.0 * (*Vabij)["acik"] * (*Tabij)["cbkj"];
  Cabij["abij"] =  2.0 * (*Vabij)["cbkj"] * (*Tabij)["acik"];
  Rabij["abij"] += Cabij["abij"];
  Rabij["abij"] += 2.0 * Cabij["acik"] * (*Tabij)["cbkj"];

  Dabij["abij"] =  (*epsi)["i"];
  Dabij["abij"] += (*epsi)["j"];
  Dabij["abij"] -= (*epsa)["a"];
  Dabij["abij"] -= (*epsa)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
  Dabij["abij"] = Dabij["abij"];

  Bivar_Function<> fDivide(&divide<double>);
  Tabij->contract(1.0, Rabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);
}

