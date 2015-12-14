#include <ParticleHoleCoulomb.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(ParticleHoleCoulomb);

ParticleHoleCoulomb::ParticleHoleCoulomb(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ParticleHoleCoulomb::~ParticleHoleCoulomb() {
}

/**
 * \brief Calculates Coulomb integrals Vabij from aiCoulombVertexReal/Imag
 */
void ParticleHoleCoulomb::run() {
  Tensor<> *aiCoulombVertexReal(
    getTensorArgument("aiCoulombVertexReal")
  );
  Tensor<> *aiCoulombVertexImag(
    getTensorArgument("aiCoulombVertexImag")
  );
  // allocate
  int nv(aiCoulombVertexReal->lens[1]);
  int no(aiCoulombVertexReal->lens[2]);
  int lens[] = { nv, nv, no, no };
  int syms[] = { NS, NS, NS, NS };
  Tensor<> *vabij(new Tensor<>(4, lens, syms, *Cc4s::world, "Cabij"));
  allocatedTensorArgument("vabij", vabij);
  (*vabij)["abij"] =  (*aiCoulombVertexReal)["gai"]*
                                 (*aiCoulombVertexReal)["gbj"];

  (*vabij)["abij"] += (*aiCoulombVertexImag)["gai"]*
                                 (*aiCoulombVertexImag)["gbj"];
// read from tensors: aiCoulombVertexImagData->value
// allocate and write to tensor vabijData->value
//  (*vabijData->value)["i.."] = (*
}

