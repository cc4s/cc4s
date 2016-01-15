#include <DrccdEnergyFromCoulombVertex.hpp>
#include <util/Log.hpp>
#include <util/MathFunctions.hpp>
#include <util/ComplexTensor.hpp>
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
 * \brief Calculates MP2 energy from Coulomb integrals Vabij
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
  Tensor<> *realTabij(new Tensor<>(4, lens, syms, *Cc4s::world, "Tabij"));
  allocatedTensorArgument("DrccdDoublesAmplitudes", realTabij);
  vabij = new Tensor<complex>(
    4, realTabij->lens, realTabij->sym, *Cc4s::world, "Vabij"
  );

  Tensor<complex> conjGammaGai(*gammaGai);
  Univar_Function<complex> fConj(&conj<complex>);
  conjGammaGai.sum(1.0, *gammaGai,"Gai", 0.0,"Gai", fConj);
  (*vabij)["abij"] = conjGammaGai["Gai"] * (*gammaGai)["Gbj"];
  Tensor<> realVabij(false, *realTabij);
  Tensor<> imagVabij(false, *realTabij);
  fromComplexTensor(*vabij, realVabij, imagVabij);
  double error(imagVabij.norm2());
  LOG(4) << "|imag(Vabij)| = " << error << std::endl;

  Scalar<> energy(*Cc4s::world);
  double e(0), dire, exce;

  LOG(0) <<
    "Solving direct ring Coupled Cluster Doubles Amplitude Equations:" <<
    std::endl;
 
  for (int i(0); i < Cc4s::options->niter; ++i) {
    LOG(0) << "iteration: " << i << std::endl;
    iterate();
    energy[""] = 2.0 * (*realTabij)["abij"] * realVabij["abij"];
    dire = energy.get_val();
    energy[""] = -1.0 * (*realTabij)["abji"] * realVabij["abij"];
    exce = energy.get_val();
    e = dire + exce;
    LOG(0) << "e=" << e << std::endl;
    LOG(1) << "RPA=" << dire << std::endl;
    LOG(1) << "SOSEX=" << exce << std::endl;
  }

  setRealArgument("DrccdEnergy", e);
}

void DrccdEnergyFromCoulombVertex::iterate() {
  // get tensors
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<complex> *gammaGai(
    getTensorArgument<complex>("ParticleHoleCoulombVertex")
  );
  Tensor<> *realTabij(getTensorArgument("DrccdDoublesAmplitudes"));
  Tensor<> imagTabij(false, *realTabij);
  Tensor<complex> tabij(
    4, realTabij->lens, realTabij->sym, *Cc4s::world, "Tabij"
  );
  toComplexTensor(*realTabij, imagTabij, tabij);
  Tensor<complex> Rabij(false, tabij);
  Tensor<complex> conjGammaGai(false, *gammaGai);
  Univar_Function<complex> fConj(&conj<complex>);
  conjGammaGai.sum(1.0, *gammaGai,"Gai", 0.0,"Gai", fConj);


  // Coulomb term
  Rabij["abij"] = (*vabij)["abij"];

  // Bubble on the right
  Tensor<complex> LGai(false, *gammaGai);
  LGai["Gai"] = 2.0 * tabij["acik"] * conjGammaGai["Gck"];
  Rabij["abij"] += LGai["Gai"] * (*gammaGai)["Gbj"];

  // Bubble on the left
  Tensor<complex> RGai(false, *gammaGai);
  RGai["Gbj"] = 2.0 * (*gammaGai)["Gck"] * tabij["cbkj"];
  Rabij["abij"] += conjGammaGai["Gai"] * RGai["Gbj"];

  // Bubbles on both sides
  Rabij["abij"] += LGai["Gai"] * RGai["Gbj"];

  Tensor<> realRabij(false, *realTabij);
  Tensor<> imagRabij(false, *realTabij);
  fromComplexTensor(Rabij, realRabij, imagRabij);
  double error(imagRabij.norm2());
  LOG(4) << "|imag(Rabij)| = " << error << std::endl;

  // to complex
  Tensor<> Dabij(false, *realTabij);
  Dabij["abij"] =  (*epsi)["i"];
  Dabij["abij"] += (*epsi)["j"];
  Dabij["abij"] -= (*epsa)["a"];
  Dabij["abij"] -= (*epsa)["b"];

  Bivar_Function<> fDivide(&divide<double>);
  realTabij->contract(1.0, realRabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);
}

