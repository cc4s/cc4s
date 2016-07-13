#include <algorithms/ReduceParticleHoleCoulombVertex.hpp>
#include <util/DryTensor.hpp>
#include <util/Exception.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ReduceParticleHoleCoulombVertex);

ReduceParticleHoleCoulombVertex::ReduceParticleHoleCoulombVertex(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ReduceParticleHoleCoulombVertex::~ReduceParticleHoleCoulombVertex() {
}

void ReduceParticleHoleCoulombVertex::run() {
  Tensor<complex> *GammaGai(
    getTensorArgument<complex>("ParticleHoleCoulombVertex")
  );
  Tensor<complex> *UGg(
    getTensorArgument<complex>("EnergyMatrixTransform")
  );
  int lens[] = { UGg->lens[1], GammaGai->lens[1], GammaGai->lens[2] };
  int syms[] = { NS, NS, NS };
  Tensor<complex> *Gammagai = new Tensor<complex>(
    3, lens, syms, *GammaGai->wrld, "Gammagai"
  );
  allocatedTensorArgument<complex>(
    "ReducedParticleHoleCoulombVertex", Gammagai
  );
  (*Gammagai)["gai"] = (*GammaGai)["Gai"] * (*UGg)["Gg"];
}

void ReduceParticleHoleCoulombVertex::dryRun() {
  DryTensor<complex> *GammaGai(
    getTensorArgument<complex, DryTensor<complex>>("ParticleHoleCoulombVertex")
  );
  DryTensor<complex> *UGg(
    getTensorArgument<complex, DryTensor<complex>>("EnergyMatrixTransform")
  );
  int lens[] = { UGg->lens[1], GammaGai->lens[1], GammaGai->lens[2] };
  int syms[] = { NS, NS, NS };
  DryTensor<complex> *Gammagai = new DryTensor<complex>(3, lens, syms);
  allocatedTensorArgument<complex, DryTensor<complex>>(
    "ReducedParticleHoleCoulombVertex", Gammagai
  );
}

