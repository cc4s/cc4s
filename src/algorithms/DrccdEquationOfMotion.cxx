#include <algorithms/DrccdEquationOfMotion.hpp>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <util/CtfMachineTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace tcc;

ALGORITHM_REGISTRAR_DEFINITION(DrccdEquationOfMotion);

DrccdEquationOfMotion::DrccdEquationOfMotion(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

DrccdEquationOfMotion::~DrccdEquationOfMotion() {
}

void DrccdEquationOfMotion::run() {
  typedef CtfMachineTensor<> MT;
  typedef CTF::Tensor<> T;

  auto machineTensorFactory(MT::Factory::create());
  auto tcc(Tcc<double>::create(machineTensorFactory));

  // Read the Drccd doubles amplitudes Tabij
  T *ctfTabij( getTensorArgument<double, T>("DrccdDoublesAmplitudes") );
  auto Tabij( tcc->createTensor(MT::create(*ctfTabij)) );

  // Read the required Coulomb integrals Vabij
  T *ctfVabij( getTensorArgument<double, T>("PPHHCoulombIntegrals") );
  auto Vabij( tcc->createTensor(MT::create(*ctfVabij)) );

  // Read and bulld energy denominators
  T *ctfEpsi( getTensorArgument<>("HoleEigenEnergies") );
  T *ctfEpsa( getTensorArgument<>("ParticleEigenEnergies") );
  // same for the HF eigenvalues
  auto epsi( tcc->createTensor(MT::create(*ctfEpsi)) );
  auto epsa( tcc->createTensor(MT::create(*ctfEpsa)) );

  // convert energy denominators tensor into tcc tensor
  auto Dabij( tcc->createTensor(Tabij, "Dabij") );

  tcc->compile( (
    (*Dabij)["abij"] <<= (*epsa)["a"],
    (*Dabij)["abij"]  += (*epsa)["b"],
    (*Dabij)["abij"]  -= (*epsi)["i"],
    (*Dabij)["abij"]  -= (*epsi)["j"]
  ) )->execute();

  // create Right eigenvector amplitudes Rabij and residuum Xabij of same shape
  // as doubles amplitudes
  auto Rabij( tcc->createTensor(Tabij, "Rabij") );
  auto Xabij( tcc->createTensor(Tabij, "Xabij") );

  // Similarity transformed hamiltonian, e^{-T} H e^{T}
  auto Habij( tcc->createTensor(Tabij, "Habij") );

  double spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);

  // Build the non-trivial part of the similarity transformed hamiltonian
  // \bar H
  tcc->compile( (
    (*Habij)["abij"] <<= (*Vabij)["abij"],
    (*Habij)["abij"]  += spins * (*Vabij)["acik"] * (*Tabij)["cbkj"]
  ) )->execute();

  // Build \bar H R
  auto iterationOperation(
    tcc->compile( (
      (*Xabij)["abij"] <<= (*Dabij)["abij"] * (*Rabij)["abij"],
      (*Xabij)["abij"]  += (*Rabij)["acik"] * (*Habij)["cbkj"],
      (*Xabij)["abij"]  += (*Habij)["caki"] * (*Rabij)["bcjk"]
    ) )
  );

  // execute iterations
  int maxIterations(getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS));
  for (int i(0); i < maxIterations; ++i) {
    LOG(0, "R_EquationOfMotion") << "Iteration " << i << "..." << std::endl;
    iterationOperation->execute();
  }

  //allocatedTensorArgument<double, T>(
    //"DrccdLambdaDoublesAmplitudes",
    //new T(Labij->template getMachineTensor<MT>()->tensor)
  //);

}
