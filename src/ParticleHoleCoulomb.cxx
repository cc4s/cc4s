#include <ParticleHoleCoulomb.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace CTF;

ParticleHoleCoulomb::ParticleHoleCoulomb(
  std::vector<Argument const *> const &argumentList
): Algorithm(argumentList) {
  
}

ParticleHoleCoulomb::~ParticleHoleCoulomb() {
}

/**
 * \brief Calculates Coulomb integrals Vabij from aiCoulombVertexReal/Imag
 */
void ParticleHoleCoulomb::run() {
  TensorData<> *aiCoulombVertexRealData(
    getTensorDataArgument("aiCoulombVertexReal")
  );
  TensorData<> *aiCoulombVertexImagData(
    getTensorDataArgument("aiCoulombVertexImag")
  );
  TensorData<> *vabijData(getTensorDataArgument("vabij"));
  // allocate
  int nv(aiCoulombVertexRealData->value->lens[1]);
  int no(aiCoulombVertexRealData->value->lens[2]);
  int lens[] = { nv, nv, no, no };
  int syms[] = { NS, NS, NS, NS };
  vabijData->value = new Tensor<>(4, lens, syms, *Cc4s::world, "Cabij");
 
// read from tensors: aiCoulombVertexImagData->value
// allocate and write to tensor vabijData->value
//  (*vabijData->value)["i.."] = (*
}

