#include <algorithms/PolarizabilityFromCoulombVertex.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <util/DryTensor.hpp>
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
 
  Matrix<> Deltaai(epsa->lens[0], epsi->lens[1], NS, *epsa->wrld, "Deltaai");
  Deltaai["ai"] =  (*epsa)["a"];
  Deltaai["ai"] -= (*epsi)["i"];

  Tensor<complex> *GammaGai(
    getTensorArgument<complex>("ParticleHoleCoulombVertex")
  );
  // calculate conjugate of GammaGai
  Tensor<complex> conjGammaGai(false, *GammaGai);
  Univar_Function<complex> fConj(&conj<complex>);
  conjGammaGai.sum(1.0, *GammaGai,"Gai", 0.0,"Gai", fConj);

  Matrix<complex> *Chi = new Matrix<complex>(
    GammaGai->lens[0], GammaGai->lens[0], *GammaGai->wrld
  );
  allocatedTensorArgument<complex>("StaticPolarizability", energyMatrix);

  (*energyMatrix)["GH"] =
    conjGammaGai["Gai"] * (*GammaGai)["Hbj"] * cTabij["abij"];
  (*energyMatrix)["GH"] *= 2.0;

  Scalar<complex> energy(*energyMatrix->wrld);
  energy[""] = (*energyMatrix)["GG"];
  complex e(energy.get_val());
  LOG(1, "EMAT") << "Tr{E}=" << e << std::endl;
}

void PolarizabilityFromCoulombVertex::dryRun() {
  DryTensor<> *Tabij(
    getTensorArgument<double, DryTensor<double>>("DoublesAmplitudes")
  );
 
  // TODO: use complex conversion routines
  DryTensor<> imagTabij(*Tabij);
  DryTensor<complex> cTabij(4, Tabij->lens.data(), Tabij->syms.data());

  DryTensor<complex> *GammaGai(
    getTensorArgument<complex, DryTensor<complex>>("ParticleHoleCoulombVertex")
  );

  // calculate conjugate of GammaGai
  DryTensor<complex> conjGammaGai(*GammaGai);

  DryMatrix<complex> *energyMatrix = new DryMatrix<complex>(
    GammaGai->lens[0], GammaGai->lens[0], NS
  );
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "EnergyMatrix", energyMatrix
  );

  // Allocate intermediate (implicitly used)
  int syms[] = { NS, NS, NS };
  int Gai[] = { GammaGai->lens[0], Tabij->lens[0], Tabij->lens[2] };
  DryTensor<complex> GammaTGai(3, Gai, syms);
}

