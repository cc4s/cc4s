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
  Tensor<> *vabij(getTensorArgument("ParticleHoleCoulombIntegrals"));

  int nv(vabij->lens[0]);
  int no(vabij->lens[2]);
  int lens[] = { nv, nv, no, no };
  int syms[] = { NS, NS, NS, NS };
  Tensor<> *tabij(new Tensor<>(4, lens, syms, *Cc4s::world, "Tabij"));
  allocatedTensorArgument("DrccdDoublesAmplitudes", tabij);
  
  Scalar<> energy(*Cc4s::world);
  double e(0), dire, exce;

  LOG(0) <<
    "Solving direct ring Coupled Cluster Doubles Amplitude Equations:" <<
    std::endl;
 
  for (int i(0); i < Cc4s::options->niter; ++i) {
    LOG(0) << "iteration: " << i << std::endl;
    iterate();
    energy[""] = 2.0 * (*tabij)["abij"] * (*vabij)["abij"];
    dire = energy.get_val();
    energy[""] = (*tabij)["abji"] * (*vabij)["abij"];
    exce = -1.0 * energy.get_val();
    e = dire + exce;
    LOG(0) << "e=" << e << std::endl;
  }

  setRealArgument("DrccdEnergy", e);
}

void DrccdEnergyFromCoulombIntegrals::iterate() {
  // get tensors
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<> *vabij(getTensorArgument("ParticleHoleCoulombIntegrals"));
  Tensor<> *tabij(getTensorArgument("DrccdDoublesAmplitudes"));
  Tensor<> Rabij(*vabij);
  Tensor<> Cabij(*vabij);
  Tensor<> Dabij(*vabij);

  Rabij["abij"] = (*vabij)["abij"];
  Rabij["abij"] += 2.0 * (*vabij)["acik"] * (*tabij)["cbkj"];
  Cabij["abij"] =  2.0 * (*vabij)["cbkj"] * (*tabij)["acik"];
  Rabij["abij"] += Cabij["abij"];
  Rabij["abij"] += 2.0 * Cabij["acik"] * (*tabij)["cbkj"];

  Dabij["abij"] += (*epsi)["i"];
  Dabij["abij"] += (*epsi)["j"];
  Dabij["abij"] -= (*epsa)["a"];
  Dabij["abij"] -= (*epsa)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
  Dabij["abij"] = Dabij["abij"];

  Bivar_Function<> fDivide(&divide<double>);
  tabij->contract(1.0, Rabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);
}

