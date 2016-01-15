#include <ParticleHoleCoulombIntegrals.hpp>
#include <util/Complex.hpp>
#include <util/ComplexTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(ParticleHoleCoulombIntegrals);

ParticleHoleCoulombIntegrals::ParticleHoleCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ParticleHoleCoulombIntegrals::~ParticleHoleCoulombIntegrals() {
}

/**
 * \brief Calculates Coulomb integrals from aiCoulombVertexReal/Imag
 */
void ParticleHoleCoulombIntegrals::run() {
  Tensor<complex> *gammaGai(
    getTensorArgument<complex>("ParticleHoleCoulombVertex")
  );
  Tensor<> realGammaGai(
    3, gammaGai->lens, gammaGai->sym, *gammaGai->wrld, "RealGammaGai"
  );
  Tensor<> imagGammaGai(
    3, gammaGai->lens, gammaGai->sym, *gammaGai->wrld, "ImagGammaGai"
  );
  // split into real and imaginary parts
  fromComplexTensor(*gammaGai, realGammaGai, imagGammaGai);

  // allocate coulomb integrals
  int nv(gammaGai->lens[1]);
  int no(gammaGai->lens[2]);
  int lens[] = { nv, nv, no, no };
  int syms[] = { NS, NS, NS, NS };
  Tensor<> *vabij(new Tensor<>(4, lens, syms, *Cc4s::world, "Vabij"));
  allocatedTensorArgument("ParticleHoleCoulombIntegrals", vabij);
  (*vabij)["abij"] =  realGammaGai["gai"] * realGammaGai["gbj"];
  (*vabij)["abij"] += imagGammaGai["gai"] * imagGammaGai["gbj"];
}

