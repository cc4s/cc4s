#include <algorithms/CoulombVertexFromFactors.hpp>
#include <math/Complex.hpp>
#include <tcc/Tcc.hpp>
#include <util/CtfMachineTensor.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

#include <vector>
#include <string>
#include <memory>

using namespace cc4s;
using namespace tcc;
using std::shared_ptr;
using std::make_shared;

ALGORITHM_REGISTRAR_DEFINITION(CoulombVertexFromFactors);

CoulombVertexFromFactors::CoulombVertexFromFactors(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

CoulombVertexFromFactors::~CoulombVertexFromFactors() {
}

void CoulombVertexFromFactors::run() {
  run<CTF::Tensor<complex>, CtfMachineTensor<complex>>(false);
}

void CoulombVertexFromFactors::dryRun() {
  run<DryTensor<complex>, DryMachineTensor<complex>>(true);
}

// TMT is either CtfMachineTensor or DryMachineTensor
template <typename T, typename MT>
void CoulombVertexFromFactors::run(const bool dryRun) {
  // create the MachineTensor factory for constructing intermediate results
  shared_ptr<MachineTensorFactory<complex>> machineTensorFactory(
    MT::Factory::create()
  );
  // create the tensor contraction compiler object
  shared_ptr<Tcc<complex>> tcc(Tcc<complex>::create(machineTensorFactory));


  // Read the Coulomb vertex GammaGqr
  T *ctfPirR( getTensorArgument<complex, T>("FactorOrbitals") );
  T *ctfLambdaFR( getTensorArgument<complex, T>("CoulombFactors") );

  // for now: create tcc::Tensors from them
  // later there will only be tcc::Tensors objects stored in cc4s
  shared_ptr<Tensor<complex>> PirR(
    tcc->createTensor(make_shared<MT>(*ctfPirR))
  );
  shared_ptr<Tensor<complex>> LambdaFR(
    tcc->createTensor(make_shared<MT>(*ctfLambdaFR))
  );

  // allocate tcc::Tensor for final result
  int NF(LambdaFR->lens[0]);
  int Np(PirR->lens[0]);
  shared_ptr<Tensor<complex>> GammaFqr(
    tcc->createTensor(std::vector<int>({NF,Np,Np}), "Gamma")
  );

  // compile and execute in one
  compile(
    (*GammaFqr)["Fqr"] <<= (*LambdaFR)["FR"] * (*PirR)["qR"] * (*PirR)["rR"]
  )->execute(
  );

  // for now: duplicate result
  // later Gamma will already be the object stored in cc4s
  shared_ptr<MT> implementationGammaFqr(
    std::dynamic_pointer_cast<MT>(
      GammaFqr->getMachineTensor()
    )
  );
  allocatedTensorArgument<complex, T>(
    "CoulombVertex", new T(implementationGammaFqr->tensor)
  );
}

