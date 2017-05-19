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

// TODO: Study the requirements to treat the problem with real numbers or
// complex

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
  // prevRabij of same shape as doubles amplitudes
  auto Rabij( tcc->createTensor(Tabij, "Rabij") );
  auto Labij( tcc->createTensor(Tabij, "Labij") );
  auto prevRabij( tcc->createTensor(Tabij, "prevRabij") );
  auto prevLabij( tcc->createTensor(Tabij, "prevLabij") );
  auto Xabij( tcc->createTensor(Tabij, "Xabij") );
  auto Yabij( tcc->createTensor(Tabij, "Yabij") );

  auto ctfXabij(&Xabij->getMachineTensor<MT>()->tensor);
  auto ctfRabij(&Rabij->getMachineTensor<MT>()->tensor);
  auto ctfLabij(&Labij->getMachineTensor<MT>()->tensor);

  // Similarity transformed hamiltonian, e^{-T} H e^{T}
  auto Habij( tcc->createTensor(Tabij, "Habij") );

  double spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);

  // Build the non-trivial part of the similarity transformed hamiltonian
  // \bar H_aj^ib (for real orbitals)
  tcc->compile( (
    (*Habij)["abij"] <<= (*Vabij)["abij"],
    (*Habij)["abij"]  += spins * (*Vabij)["acik"] * (*Tabij)["cbkj"]
  ) )->execute();

  determineEnergyShift();

  auto beta(tcc->createTensor(std::vector<int>(), "beta"));
  auto delta(tcc->createTensor(std::vector<int>(), "delta"));
  auto energy(tcc->createTensor(std::vector<int>(), "energy"));
  CTF::Scalar<> ctfBeta(*Cc4s::world);
  CTF::Scalar<> ctfDelta(*Cc4s::world);
  CTF::Scalar<> ctfEnergy(*Cc4s::world);

  auto buildHbarR(
    tcc->compile( (
      (*Xabij)["abij"] <<= (*Dabij)["abij"] * (*Rabij)["abij"],
      (*Xabij)["abij"]  += spins * (*Rabij)["acik"] * (*Habij)["cbkj"],
      (*Xabij)["abij"]  += spins * (*Habij)["caki"] * (*Rabij)["bcjk"],
      (*Xabij)["abij"]  -= energyShift * (*Rabij)["abij"]
    ) )
  );

  auto buildLHbarR(
    tcc->compile( (
      (*energy)[""] <<= 0.5*spins*spins * (*Xabij)["abij"] * (*Labij)["abij"],
      (*energy)[""]  -= 0.5 * spins * (*Xabij)["abij"] * (*Labij)["abji"]
    ) )
  );

  auto buildNewR(
    tcc->compile( (
      (*Yabij)["abij"] <<= -1.0*(*beta)[""] * (*prevRabij)["abij"],
      (*Yabij)["abij"]  += (*Xabij)["abij"],
      (*Yabij)["abij"]  -= (*energy)[""] * (*Rabij)["abij"],
      //TODO: Swap pointers instead of values
      (*prevRabij)["abij"] <<= (*Rabij)["abij"],
      (*Rabij)["abij"] <<= (*Yabij)["abij"]
    ) )
  );

  auto buildLHbar(
    tcc->compile( (
      (*Xabij)["abij"] <<= (*Vabij)["abij"],
      (*Xabij)["abij"]  += (*Dabij)["abij"] * (*Labij)["abij"],
      (*Xabij)["abij"]  += spins * (*Vabij)["acik"] * (*Habij)["bcjk"],
      (*Xabij)["abij"]  += spins * (*Habij)["acik"] * (*Labij)["bcjk"],
      (*Xabij)["abij"]  -= energyShift * (*Labij)["abij"]
    ) )
  );

  auto buildNewL(
    tcc->compile( (
      (*Yabij)["abij"] <<= -1.0*(*delta)[""] * (*prevLabij)["abij"],
      (*Yabij)["abij"]  += (*Xabij)["abij"],
      (*Yabij)["abij"]  -= (*energy)[""] * (*Labij)["abij"],
      //TODO: Swap pointers instead of values
      (*prevLabij)["abij"] <<= (*Labij)["abij"],
      (*Labij)["abij"] <<= (*Yabij)["abij"]
    ) )
  );

  auto buildBeta(
    tcc->compile( (
      (*beta)[""] <<= 0.5*spins*spins*(*Labij)["abij"] * (*Rabij)["abij"],
      (*beta)[""]  -= 0.5*spins*(*Labij)["abij"] * (*Rabij)["abji"]
    ) )
  );

  setRandomTensor(*ctfRabij);
  (*ctfLabij)["abij"] = (*ctfRabij)["abij"];
  buildBeta->execute();
  ctfBeta[""] = beta->getMachineTensor<MT>()->tensor[""];
  (*ctfLabij)["abij"] = ( 1 / ctfBeta.get_val() ) * (*ctfLabij)["abij"];

  // execute iterations
  int maxIterations(getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS));
  for (int i(0); i < maxIterations; ++i) {
    LOG(0, "DrccdEOM") << "Iteration " << i << "..." << std::endl;

    buildHbarR->execute();
    buildLHbarR->execute();
    buildNewR->execute();
    buildLHbar->execute();
    buildNewL->execute();

    // Build coefficients
    buildBeta->execute();
    ctfBeta[""] = beta->getMachineTensor<MT>()->tensor[""];
    ctfDelta[""] = std::sqrt(std::abs( ctfBeta.get_val() ));
    ctfBeta[""] = ctfBeta.get_val() / ctfDelta.get_val();

    (*ctfRabij)["abij"] = ( 1 / ctfBeta.get_val() ) * (*ctfRabij)["abij"];
    (*ctfLabij)["abij"] = ( 1 / ctfDelta.get_val() ) * (*ctfLabij)["abij"];
    beta->getMachineTensor<MT>()->tensor[""] = ctfBeta[""];
    delta->getMachineTensor<MT>()->tensor[""] = ctfDelta[""];
    double b(ctfBeta.get_val()), d(ctfDelta.get_val());
    LOG(1, "DrccdEOM") << "beta=" << b << ", delta=" << d << std::endl;

    ctfEnergy[""] = energy->getMachineTensor<MT>()->tensor[""];
    double e(ctfEnergy.get_val());

    LOG(0, "DrccdEOM") << "e= " << e+energyShift << std::endl;
  }

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
