#include <algorithms/DrccdEquationOfMotion.hpp>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <math/MathFunctions.hpp>
#include <math/RandomTensor.hpp>
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

  // create Right and Left eigenvector amplitudes Rabij, Labij and itermediate
  // Xabij of same shape as doubles amplitudes
  auto Rabij( tcc->createTensor(Tabij, "Rabij") );
  auto Labij( tcc->createTensor(Tabij, "Labij") );
  auto Xabij( tcc->createTensor(Tabij, "Xabij") );

  auto ctfRabij(&Rabij->getMachineTensor<MT>()->tensor);
  auto ctfLabij(&Labij->getMachineTensor<MT>()->tensor);
  auto ctfXabij(&Xabij->getMachineTensor<MT>()->tensor);

  setRandomTensor(*ctfRabij);
  //setRandomTensor(*ctfLabij);

  // Similarity transformed hamiltonian, e^{-T} H e^{T}
  auto Habij( tcc->createTensor(Tabij, "Habij") );

  double spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);

  // Build the non-trivial part of the similarity transformed hamiltonian
  // \bar H
  tcc->compile( (
    (*Habij)["abij"] <<= (*Vabij)["abij"],
    (*Habij)["abij"]  += spins * (*Vabij)["acik"] * (*Tabij)["cbkj"]
  ) )->execute();

  determineEnergyShift();

  // Build \bar H R
  auto rightIterationOperation(
    tcc->compile( (
      (*Xabij)["abij"] <<= (*Dabij)["abij"] * (*Rabij)["abij"],
      (*Xabij)["abij"]  += spins * (*Rabij)["acik"] * (*Habij)["cbkj"],
      (*Xabij)["abij"]  += spins * (*Habij)["caki"] * (*Rabij)["bcjk"],
      (*Xabij)["abij"]  -= energyShift * (*Rabij)["abij"]
    ) )
  );

  // Build L \bar H
  auto leftIterationOperation(
    tcc->compile( (
      (*Xabij)["abij"] <<= (*Vabij)["abij"],
      (*Xabij)["abij"]  += (*Dabij)["abij"] * (*Labij)["abij"],
      (*Xabij)["abij"]  += spins * (*Vabij)["acik"] * (*Habij)["cbkj"],
      (*Xabij)["abij"]  += spins * (*Habij)["caki"] * (*Labij)["bcjk"],
      (*Xabij)["abij"]  -= energyShift * (*Labij)["abij"]
    ) )
  );

  // execute iterations
  int maxIterations(getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS));
  for (int i(0); i < maxIterations; ++i) {
    LOG(0, "Right DrccdEOM") << "Iteration " << i << "..." << std::endl;

    rightIterationOperation->execute();
    double rNorm(frobeniusNorm(*ctfXabij));
    LOG(1, "DrccdEOM") << "|H'R|-Delta=" << rNorm-energyShift << std::endl;
    tcc->compile(
      (*Rabij)["abij"] <<= (1/rNorm) * (*Xabij)["abij"]
    )->execute();
  }

  (*ctfLabij)["abij"] = (*ctfRabij)["abij"];

  for (int i(0); i < maxIterations; ++i) {
    LOG(0, "Left DrccdEOM") << "Iteration " << i << "..." << std::endl;

    leftIterationOperation->execute();
    double lNorm(frobeniusNorm(*ctfXabij));
    LOG(1, "DrccdEOM") << "|LH'|-Delta=" << lNorm-energyShift << std::endl;
    tcc->compile(
      (*Labij)["abij"] <<= (1/lNorm) * (*Xabij)["abij"]
    )->execute();
  }

  CTF::Scalar<> biNorm;
  biNorm[""] = (*ctfRabij)["abij"] * (*ctfLabij)["abij"];

  double biNormValue(biNorm.get_val());
  LOG(1, "DrccdEOM") << "|LR| " << biNormValue << std::endl;

  (*ctfRabij)["abij"] = ( 1 / std::sqrt(biNormValue) ) * (*ctfRabij)["abij"];
  (*ctfLabij)["abij"] = ( 1 / std::sqrt(biNormValue) ) * (*ctfLabij)["abij"];

  // Compute Xabij with the binormalized L
  leftIterationOperation->execute();

  (*ctfXabij)["abij"] += energyShift * (*ctfLabij)["abij"];

  CTF::Scalar<> energy;

  energy[""] = 0.5 * spins * spins * (*ctfXabij)["abij"] * (*ctfRabij)["abij"];
  energy[""] -= 0.5 * spins * (*ctfXabij)["abij"] * (*ctfRabij)["abji"];

  //double e(ctfEnergy->get_val());
  LOG(0, "DrccdEOM") << "e= " << energy << std::endl;

}


void DrccdEquationOfMotion::determineEnergyShift() {
  typedef CTF::Tensor<> T;

  // Read energies
  T *epsi( getTensorArgument<>("HoleEigenEnergies") );
  T *epsa( getTensorArgument<>("ParticleEigenEnergies") );

  int64_t indicesCount;
  double *energies;

  // Get all data from epsi
  epsi->read_all(&indicesCount, &energies);
  energyShift = -energies[0];
  free(energies);

  // Get all data from epsa
  epsa->read_all(&indicesCount, &energies);
  energyShift += energies[indicesCount-1];
  free(energies);

  // many-body system can have No times single body excitation energies
  energyShift *= epsi->lens[0];

  // allow manual override
  energyShift = getRealArgument("energyShift", energyShift);

  LOG(1, "DrccdEOM") << "H'=H-Delta with Delta=" << energyShift << std::endl;
}
