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

  auto ctfRabij(&Rabij->getMachineTensor<MT>()->tensor);
  auto ctfLabij(&Labij->getMachineTensor<MT>()->tensor);

  auto ctfPrevRabij(&prevRabij->getMachineTensor<MT>()->tensor);
  auto ctfPrevLabij(&prevLabij->getMachineTensor<MT>()->tensor);

  auto ctfXabij(&Xabij->getMachineTensor<MT>()->tensor);

  setRandomTensor(*ctfRabij);

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

  auto ctfBeta(
    std::dynamic_pointer_cast<T>(std::make_shared<CTF::Scalar<>>())
  );
  auto ctfDelta(
    std::dynamic_pointer_cast<T>(std::make_shared<CTF::Scalar<>>())
  );
  auto ctfEnergy(
    std::dynamic_pointer_cast<T>(std::make_shared<CTF::Scalar<>>())
  );

  auto beta(tcc->createTensor(ctfBeta));
  auto delta(tcc->createTensor(ctfDelta));
  auto energy(tcc->createTensor(ctfEnergy));

  auto buildHbarR(
    tcc->compile( (
      (*Xabij)["abij"] <<= (*Dabij)["abij"] * (*Rabij)["abij"],
      (*Xabij)["abij"]  += spins * (*Rabij)["acik"] * (*Habij)["cbkj"],
      (*Xabij)["abij"]  += spins * (*Habij)["caki"] * (*Rabij)["bcjk"],
      (*Xabij)["abij"]  -= energyShift * (*Rabij)["abij"]
    ) )
  );

  auto buildLHbarR(
    tcc->complile( (
      (*energy)[""] <<= 0.5*spins*spins * (*Xabij)["abij"] * (*Labij)["abij"],
      (*energy)[""]  -= 0.5 * spins * (*Xabij)["abij"] * (*Labij)["abji"]
    ) )
  );

  auto buildNewR(
    tcc->compile( (
      (*prevRabij)["abij"] <<= -1*(*beta)[""] * (*prevRabij)["abij"],
      (*prevRabij)["abij"]  += (*Xabij)["abij"],
      (*prevRabij)["abij"]  -= (*alpha)[""] * (*Rabij)["abij"],
      // TODO: Swap pointers instead of values
      (*Xabij)["abij"] <<= (*prevRabij)["abij"],
      (*prevRabij)["abij"] <<= (*Rabij)["abij"],
      (*Rabij)["abij"] <<= (*Xabij)["abij"]
    ) )
  );

  auto buildLHbar(
    tcc->compile( (
      (*Xabij)["abij"] <<= (*Vabij)["abij"],
      (*Xabij)["abij"]  += (*Dabij)["abij"] * (*Labij)["abij"],
      (*Xabij)["abij"]  += spins * (*Vabij)["acik"] * (*Habij)["cbkj"],
      (*Xabij)["abij"]  += spins * (*Habij)["caki"] * (*Labij)["bcjk"],
      (*Xabij)["abij"]  -= energyShift * (*Labij)["abij"]
    ) )
  );

  auto buildNewL(
    tcc->compile( (
      (*prevLabij)["abij"] <<= -1*(*beta)[""] * (*prevLabij)["abij"],
      (*prevLabij)["abij"]  += (*Xabij)["abij"],
      (*prevLabij)["abij"]  -= (*alpha)[""] * (*Labij)["abij"],
      // TODO: Swap pointers instead of values
      (*Xabij)["abij"] <<= (*prevLabij)["abij"],
      (*prevLabij)["abij"] <<= (*Labij)["abij"],
      (*Labij)["abij"] <<= (*Xabij)["abij"]
    ) )
  );

  auto buildBeta(
    tcc->compile( (
      (*beta)[""] <<= 0.5*spins*spins*(*Labij)["abij"] * (*Rabij)["abij"],
      (*beta)[""]  -= 0.5*spins*(*Labij)["abij"] * (*Rabij)["abji"]
    ) )
  );

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
    ctfDelta->set_val(std::sqrt(std::abs(ctfBeta->get_val())));
    ctfBeta->set_val(
      ctfBeta->get_val() / ctfDelta->get_val()
    );

    ctfRabij["abij"] = ( 1 / ctfBeta->get_val() ) * ctfRabij["abij"];
    ctfLabij["abij"] = ( 1 / ctfDelta->get_val() ) * ctfLabij["abij"];

    double e(ctfEnergy->get_val());

    LOG(0, "DrccdEOM") << "e= " << e << std::endl;

  }

  //CTF::Scalar<> biNorm;
  //biNorm[""] = (*ctfRabij)["abij"] * (*ctfLabij)["abij"];

  //double biNormValue(biNorm.get_val());
  //LOG(1, "DrccdEOM") << "|LR| " << biNormValue << std::endl;

  //double lNorm(ctfLabij->norm1());
  //double rNorm(ctfRabij->norm1());
  //LOG(1, "DrccdEOM") << "|L| " << lNorm << std::endl;
  //LOG(1, "DrccdEOM") << "|R| " << rNorm << std::endl;

  //(*ctfRabij)["abij"] = ( 1 / std::sqrt(biNormValue) ) * (*ctfRabij)["abij"];
  //(*ctfLabij)["abij"] = ( 1 / std::sqrt(biNormValue) ) * (*ctfLabij)["abij"];

  //// Compute prevRabij with the binormalized L
  //LHbar->execute();

  //(*ctfPrevRabij)["abij"] += energyShift * (*ctfLabij)["abij"];

  //LOG(0, "DrccdEOM") << "e= " << energy << std::endl;

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
