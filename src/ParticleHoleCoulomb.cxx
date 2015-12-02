#include <ParticleHoleCoulomb.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace cc4s;

ParticleHoleCoulomb::ParticleHoleCoulomb(
  std::vector<Argument const *> const &argumentList
): Algorithm(argumentList) {
  
}

ParticleHoleCoulombVertexMp2::~ParticleHoleCoulombVertexMp2() {
}

/**
 * \brief Calculates Coulomb integrals Vabij from aiCoulombVertexReal/Imag
 */
void ParticleHoleCoulombVertexMp2::run() {
  TensorData<> *aiCoulombVertexRealData(
    getTensorDataArgument("aiCoulombVertexReal")
  );
  TensorData<> *aiCoulombVertexImagData(
    getTensorDataArgument("aiCoulombVertexImag")
  );
  TensorData<> *vabijData(getTensorDataArgument("vabij"));
  // allocate
 Tensor<>*vabijData->value(V->abij);
 
// read from tensors: aiCoulombVertexImagData->value
// allocate and write to tensor vabijData->value
}

