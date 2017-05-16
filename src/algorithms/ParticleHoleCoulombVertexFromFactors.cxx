#include <algorithms/ParticleHoleCoulombVertexFromFactors.hpp>
#include <math/Complex.hpp>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <util/CtfMachineTensor.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

#include <vector>
#include <string>
#include <memory>

using namespace cc4s;
using namespace tcc;
using std::make_shared;

ALGORITHM_REGISTRAR_DEFINITION(ParticleHoleCoulombVertexFromFactors);

ParticleHoleCoulombVertexFromFactors::ParticleHoleCoulombVertexFromFactors(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ParticleHoleCoulombVertexFromFactors::~ParticleHoleCoulombVertexFromFactors() {
}

void ParticleHoleCoulombVertexFromFactors::run() {
  run<CTF::Tensor<complex>, CtfMachineTensor<complex>>(false);
}

void ParticleHoleCoulombVertexFromFactors::dryRun() {
  run<DryTensor<complex>, DryMachineTensor<complex>>(true);
}

// TMT is either CtfMachineTensor or DryMachineTensor
template <typename T, typename MT>
void ParticleHoleCoulombVertexFromFactors::run(const bool dryRun) {
  auto machineTensorFactory(MT::Factory::create());
  auto tcc(Tcc<complex>::create(machineTensorFactory));


  // Read the Coulomb vertex GammaGqr
  T *ctfPiiR( getTensorArgument<complex, T>("HoleFactorOrbitals") );
  T *ctfPiaR( getTensorArgument<complex, T>("ParticleFactorOrbitals") );
  T *ctfLambdaFR( getTensorArgument<complex, T>("CoulombFactors") );

  // for now: create tcc::Tensors from them
  // later there will only be tcc::Tensors objects stored in cc4s
  auto PiiR( tcc->createTensor(MT::create(*ctfPiiR)) );
  auto PiaR( tcc->createTensor(MT::create(*ctfPiaR)) );
  auto LambdaFR( tcc->createTensor(MT::create(*ctfLambdaFR)) );
 
  // allocate tcc::Tensor for final result
  int NF(LambdaFR->lens[0]);
  int Nv(PiaR->lens[0]);
  int No(PiiR->lens[0]);
  auto GammaFai( tcc->createTensor(std::vector<int>({NF,Nv,No}), "GammaFai") );

  // compile
  auto operation(
    tcc->compile(
      (*GammaFai)["Fai"] <<= (*LambdaFR)["FR"] * (*PiaR)["aR"] * (*PiiR)["iR"]
    )
  );
  // and execute
  operation->execute();

  // for now: duplicate result
  // later Gamma will already be the object stored in cc4s
  allocatedTensorArgument<complex, T>(
    "ParticleHoleCoulombVertex",
    new T(GammaFai->template getMachineTensor<MT>()->tensor)
  );
}

