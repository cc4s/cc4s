#include <DrccdEnergyFromCoulombVertex.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <Options.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DrccdEnergyFromCoulombVertex);

DrccdEnergyFromCoulombVertex::DrccdEnergyFromCoulombVertex(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

DrccdEnergyFromCoulombVertex::~DrccdEnergyFromCoulombVertex() {
}

/**
 * \brief Calculates Drccd energy from Coulomb vertex gammaGai
 */
void DrccdEnergyFromCoulombVertex::run() {
  Tensor<complex> *gammaGai(
    getTensorArgument<complex>("ParticleHoleCoulombVertex")
  );
  realGammaGai = new Tensor<>(
    3, gammaGai->lens, gammaGai->sym, *Cc4s::world, "realGammaGai"
  );
  imagGammaGai = new Tensor<>(
    3, gammaGai->lens, gammaGai->sym, *Cc4s::world, "imagGammaGai"
  );
  fromComplexTensor(*gammaGai, *realGammaGai, *imagGammaGai);

  int nv(gammaGai->lens[1]);
  int no(gammaGai->lens[2]);
  int lens[] = { nv, nv, no, no };
  int syms[] = { NS, NS, NS, NS };
  Tensor<> *tabij(new Tensor<>(4, lens, syms, *Cc4s::world, "Tabij"));
  allocatedTensorArgument("DrccdDoublesAmplitudes", tabij);
  vabij = new Tensor<>(4, tabij->lens, tabij->sym, *Cc4s::world, "Vabij");

  (*vabij)["abij"] =  (*realGammaGai)["Gai"] * (*realGammaGai)["Gbj"];
  (*vabij)["abij"] += (*imagGammaGai)["Gai"] * (*imagGammaGai)["Gbj"];

  // NOTE: for debugging:
  Tensor<> imagVabij(false, *tabij);
  imagVabij["abij"] =  (*realGammaGai)["Gai"] * (*imagGammaGai)["Gbj"];
  imagVabij["abij"] -= (*imagGammaGai)["Gai"] * (*realGammaGai)["Gbj"];
  double error(imagVabij.norm2());
  LOG(4) << "|imag(Vabij)| = " << error << std::endl;

  // allocate intermediate tensors
  Rabij = new Tensor<>(false, *tabij);
  Dabij = new Tensor<>(false, *tabij);
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  (*Dabij)["abij"] =  (*epsi)["i"];
  (*Dabij)["abij"] += (*epsi)["j"];
  (*Dabij)["abij"] -= (*epsa)["a"];
  (*Dabij)["abij"] -= (*epsa)["b"];
  realLGai = new Tensor<>(false, *realGammaGai);
  imagLGai = new Tensor<>(false, *imagGammaGai);
  realRGai = new Tensor<>(false, *realGammaGai);
  imagRGai = new Tensor<>(false, *imagGammaGai);

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
    energy[""] = -1.0 * (*tabij)["abji"] * (*vabij)["abij"];
    exce = energy.get_val();
    e = dire + exce;
    LOG(0) << "e=" << e << std::endl;
    LOG(1) << "RPA=" << dire << std::endl;
    LOG(1) << "SOSEX=" << exce << std::endl;
  }

  setRealArgument("DrccdEnergy", e);

  // deallocate intermediate tensors
  delete Rabij; delete Dabij;
  delete realLGai; delete imagLGai;
  delete realRGai; delete imagRGai;
  delete realGammaGai; delete imagGammaGai;
}

void DrccdEnergyFromCoulombVertex::iterate() {
  // get tensors
  Tensor<> *tabij(getTensorArgument("DrccdDoublesAmplitudes"));

  // Coulomb term
  (*Rabij)["abij"] = (*vabij)["abij"];

  // Bubble on the right
  // LGai["Gai"] = 2.0 * tabij["acik"] * conj(gammaGai["Gck"]);
  (*realLGai)["Gai"] = +2.0 * (*tabij)["acik"] * (*realGammaGai)["Gck"];
  (*imagLGai)["Gai"] = -2.0 * (*tabij)["acik"] * (*imagGammaGai)["Gck"];
  // Rabij["abij"] += real( LGai["Gai"] * (*gammaGai)["Gbj"] );
  (*Rabij)["abij"] += (*realLGai)["Gai"] * (*realGammaGai)["Gbj"];
  (*Rabij)["abij"] -= (*imagLGai)["Gai"] * (*imagGammaGai)["Gbj"];

  // Bubble on the left
  // RGai["Gbj"] = 2.0 * (*gammaGai)["Gck"] * tabij["cbkj"];
  (*realRGai)["Gbj"] = +2.0 * (*realGammaGai)["Gck"] * (*tabij)["cbkj"];
  (*imagRGai)["Gbj"] = +2.0 * (*imagGammaGai)["Gck"] * (*tabij)["cbkj"];
  // Rabij["abij"] += real( conj(gammaGai["Gai"]) * RGai["Gbj"] );
  (*Rabij)["abij"] += (*realGammaGai)["Gai"] * (*realRGai)["Gbj"];
  (*Rabij)["abij"] += (*imagGammaGai)["Gai"] * (*imagRGai)["Gbj"];

  // Bubbles on both sides
  // Rabij["abij"] += real( LGai["Gai"] * RGai["Gbj"] );
  (*Rabij)["abij"] += (*realLGai)["Gai"] * (*realRGai)["Gbj"];
  (*Rabij)["abij"] -= (*imagLGai)["Gai"] * (*imagRGai)["Gbj"];

  // This is what we fixed last time...
  Bivar_Function<> fDivide(&divide<double>);
  tabij->contract(1.0, *Rabij,"abij", *Dabij,"abij", 0.0,"abij", fDivide);
}

