#include <algorithms/ApproximateParticleHoleCoulombVertex.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ApproximateParticleHoleCoulombVertex);

ApproximateParticleHoleCoulombVertex::ApproximateParticleHoleCoulombVertex(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ApproximateParticleHoleCoulombVertex::~ApproximateParticleHoleCoulombVertex() {
}

void ApproximateParticleHoleCoulombVertex::run() {
  Tensor<complex> *GammaGai(
    getTensorArgument<complex>("FullParticleHoleCoulombVertex")
  );
  Tensor<complex> *UGF(
    getTensorArgument<complex>("ParticleHoleCoulombVertexSingularVectors")
  );
  int lens[] = { UGF->lens[1], GammaGai->lens[1], GammaGai->lens[2] };
  int syms[] = { NS, NS, NS };
  Tensor<complex> *GammaFai = new Tensor<complex>(
    3, lens, syms, *GammaGai->wrld, "GammaFai"
  );
  allocatedTensorArgument<complex>(
    "ParticleHoleCoulombVertex", GammaFai
  );
  (*GammaFai)["Fai"] = (*GammaGai)["Gai"] * (*UGF)["GF"];
}

void ApproximateParticleHoleCoulombVertex::dryRun() {
  DryTensor<complex> *GammaGai(
    getTensorArgument<complex, DryTensor<complex>>("FullParticleHoleCoulombVertex")
  );
  DryTensor<complex> *UGF(
    getTensorArgument<complex, DryTensor<complex>>(
      "ParticleHoleCoulombVertexSingularVectors"
    )
  );
  int lens[] = { UGF->lens[1], GammaGai->lens[1], GammaGai->lens[2] };
  int syms[] = { NS, NS, NS };
  DryTensor<complex> *GammaFai = new DryTensor<complex>(
    3, lens, syms, SOURCE_LOCATION
  );
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "ParticleHoleCoulombVertex", GammaFai
  );
}

