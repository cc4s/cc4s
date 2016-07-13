#include <algorithms/SliceCoulombVertex.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <util/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(SliceCoulombVertex);

SliceCoulombVertex::SliceCoulombVertex(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

SliceCoulombVertex::~SliceCoulombVertex() {
}

void SliceCoulombVertex::run() {
  // Read the Coulomb vertex GammaGpq
  Tensor<complex> *GammaGpq( getTensorArgument<complex>("CoulombVertex"));

  // Read the Particle/Hole Eigenenergies
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));

  // Compute the No,Nv,NG,Np
  int NG(GammaGpq->lens[0]);
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np = No + Nv;

  // Allocate and compute GammaGai
  int GaiStart[] = {0 ,No, 0};
  int GaiEnd[]   = {NG,Np,No};
  Tensor<complex> *GammaGai(
    new Tensor<complex>(GammaGpq->slice(GaiStart,GaiEnd))
  );
  allocatedTensorArgument<complex>("ParticleHoleCoulombVertex", GammaGai);
}

void SliceCoulombVertex::dryRun() {
  // Read the Coulomb vertex GammaGpq
  DryTensor<complex> *GammaGpq(
    getTensorArgument<complex, DryTensor<complex>>("CoulombVertex")
  );

  // Read the Particle/Hole Eigenenergies
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies")
  );

  // Compute the No,Nv,NG
  int NG(GammaGpq->lens[0]);
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Allocate and compute GammaGab,GammaGai,GammaGij from GammaGpq
  int GaiLens[] = { NG,Nv,No };
  int GaiSyms[] = { NS,NS,NS };
  DryTensor<complex> *GammaGai(new DryTensor<complex>(3, GaiLens, GaiSyms));

  allocatedTensorArgument<complex, DryTensor<complex>>(
    "ParticleHoleCoulombVertex", GammaGai
  );
}
