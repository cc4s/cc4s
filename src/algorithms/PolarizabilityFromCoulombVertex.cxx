#include <algorithms/PolarizabilityFromCoulombVertex.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(PolarizabilityFromCoulombVertex);

PolarizabilityFromCoulombVertex::PolarizabilityFromCoulombVertex(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

PolarizabilityFromCoulombVertex::~PolarizabilityFromCoulombVertex() {
}

void PolarizabilityFromCoulombVertex::run() {
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
 
  Matrix<> Deltaai(epsa->lens[0], epsi->lens[0], NS, *epsa->wrld, "Deltaai");
  Deltaai["ai"] =  (*epsa)["a"];
  Deltaai["ai"] -= (*epsi)["i"];

  Matrix<complex> cDeltaai(
    epsa->lens[0], epsi->lens[0], NS, *epsa->wrld, "cDeltaai"
  );
  toComplexTensor(Deltaai, cDeltaai);

  Tensor<complex> *GammaGai(
    getTensorArgument<complex>("ParticleHoleCoulombVertex")
  );
  // calculate conjugate of GammaGai
  Tensor<complex> conjGammaGai(false, *GammaGai);
  Univar_Function<complex> fConj(&conj<complex>);
  conjGammaGai.sum(1.0, *GammaGai,"Gai", 0.0,"Gai", fConj);
  Bivar_Function<complex> fDivide(&divide<complex>);
  conjGammaGai.contract(1.0, conjGammaGai,"Gai", cDeltaai,"ai", 0.0,"Gai", fDivide);

  Matrix<complex> *Chi = new Matrix<complex>(
    GammaGai->lens[0], GammaGai->lens[0], *GammaGai->wrld
  );
  allocatedTensorArgument<complex>("StaticPolarizability", Chi);

  (*Chi)["GH"] = conjGammaGai["Gai"] * (*GammaGai)["Hai"];
  // Chi = Chi + conj(Chi)
  Chi->sum(1.0, *Chi,"GH", 1.0,"GH", fConj);

  LOG(1, "StaticPolarizability") <<
    "NOTE: actually calculating v^1/2_G X0_G^H v^1/2^H" << std::endl;
}

void PolarizabilityFromCoulombVertex::dryRun() {
  //  DryTensor<> *epsa(
  getTensorArgument<double, DryTensor<>>("ParticleEigenEnergies");
  //  );
  //  DryTensor<> *epsi(
  getTensorArgument<double, DryTensor<>>("HoleEigenEnergies");
  //  );
 
  DryTensor<complex> *GammaGai(
    getTensorArgument<complex, DryTensor<complex>>("ParticleHoleCoulombVertex")
  );

  // calculate conjugate of GammaGai
  DryTensor<complex> conjGammaGai(*GammaGai);

  DryMatrix<complex> *Chi = new DryMatrix<complex>(
    GammaGai->lens[0], GammaGai->lens[0], NS
  );
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "StaticPolarizability", Chi
  );
}

