#include <algorithms/ParticleHoleCoulombIntegrals.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
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
  Tensor<complex> *GammaGai(
    getTensorArgument<complex>("ParticleHoleCoulombVertex")
  );
  Tensor<> realGammaGai(
    3, GammaGai->lens, GammaGai->sym, *GammaGai->wrld, "RealGammaGai"
  );
  Tensor<> imagGammaGai(
    3, GammaGai->lens, GammaGai->sym, *GammaGai->wrld, "ImagGammaGai"
  );
  // split into real and imaginary parts
  fromComplexTensor(*GammaGai, realGammaGai, imagGammaGai);

  // allocate coulomb integrals
  int Nv(GammaGai->lens[1]);
  int No(GammaGai->lens[2]);
  int lens[] = { Nv, Nv, No, No };
  int syms[] = { NS, NS, NS, NS };
  Tensor<> *Vabij(new Tensor<>(4, lens, syms, *Cc4s::world, "Vabij"));
  allocatedTensorArgument("PPHHCoulombIntegrals", Vabij);
  (*Vabij)["abij"] =  realGammaGai["gai"] * realGammaGai["gbj"];
  (*Vabij)["abij"] += imagGammaGai["gai"] * imagGammaGai["gbj"];
}

