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
  typedef CTF::Tensor<real> T;
  typedef tcc::Tcc<CtfEngine> TCC;

  // Read the Drccd doubles amplitudes Tabij
  T *ctfTabij( getTensorArgument<double, T>("DrccdDoublesAmplitudes") );
  auto Tabij( tcc::Tensor<complex,CtfEngine>::create(*ctfTabij) );

  // Read the required Coulomb integrals Vabij
  T *ctfVabij( getTensorArgument<double, T>("PPHHCoulombIntegrals") );
  auto Vabij( tcc::Tensor<complex,CtfEngine>::create(*ctfVabij) );

  // Read and bulld energy denominators
  T *ctfEpsi( getTensorArgument<>("HoleEigenEnergies") );
  T *ctfEpsa( getTensorArgument<>("ParticleEigenEnergies") );
  // same for the HF eigenvalues
  auto epsi( tcc::Tensor<complex,CtfEngine>::create(*ctfEpsi) );
  auto epsa( tcc::Tensor<complex,CtfEngine>::create(*ctfEpsa) );

  // convert energy denominators tensor into tcc tensor
  auto Dabij( TCC::tensor(Tabij, "Dabij") );
  auto ctfDabij(&Dabij->getMachineTensor()->tensor);

  IndexCounts indexCounts;
  (
    (*Dabij)["abij"] <<= (*epsa)["a"],
    (*Dabij)["abij"]  += (*epsa)["b"],
    (*Dabij)["abij"]  -= (*epsi)["i"],
    (*Dabij)["abij"]  -= (*epsi)["j"]
  )->compile(indexCounts)->execute();

  // create Right and Left eigenvector amplitudes Rabij, Labij and itermediate
  // prevRabij of same shape as doubles amplitudes
  auto Rabij( TCC::tensor(Tabij, "Rabij") );
  auto Labij( TCC::tensor(Tabij, "Labij") );
  auto prevRabij( TCC::tensor(Tabij, "prevRabij") );
  auto prevLabij( TCC::tensor(Tabij, "prevLabij") );
  auto Xabij( TCC::tensor(Tabij, "Xabij") );
  auto Yabij( TCC::tensor(Tabij, "Yabij") );

  auto ctfPrevLabij(&prevLabij->getMachineTensor()->tensor);
  auto ctfPrevRabij(&prevRabij->getMachineTensor()->tensor);
  auto ctfYabij(&Yabij->getMachineTensor()->tensor);
  auto ctfXabij(&Xabij->getMachineTensor()->tensor);
  auto ctfRabij(&Rabij->getMachineTensor()->tensor);
  auto ctfLabij(&Labij->getMachineTensor()->tensor);

  // Similarity transformed hamiltonian, e^{-T} H e^{T}
  auto Habij( TCC::tensor(Tabij, "Habij") );
  auto ctfHabij(&Habij->getMachineTensor()->tensor);

  double spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);

  // Build the non-trivial part of the similarity transformed hamiltonian
  // \bar H_aj^ib (for real orbitals)
  (
    (*Habij)["abij"] <<= (*Vabij)["abij"],
    (*Habij)["abij"]  += spins * (*Vabij)["acik"] * (*Tabij)["cbkj"]
  )->compile(indexCounts)->execute();

  if (getIntegerArgument("fullDiagonalization", 0)) {
    auto H20(new CTF::Tensor<>(*ctfVabij));

    (*H20)["abij"] -= (*ctfVabij)["abji"];

    int Nv(epsa->lens[0]), No(epsi->lens[0]);
    // cdkl are row indices and abij are column indices
    int lens[] = {Nv,Nv,No,No, Nv,Nv,No,No};
    int syms[] = {NS,NS,NS,NS, NS,NS,NS,NS};
    auto H22( new CTF::Tensor<>(8, lens, syms, *Cc4s::world, "H22cdklabij") );
    // diagonal elements
    (*H22)["abijabij"] = (*ctfDabij)["abij"];
    (*H22)["cbkjabij"] += (*ctfHabij)["acik"];
    (*H22)["adilabij"] += (*ctfHabij)["bdjl"];

    //(*H22)["abjiabij"] -= (*ctfDabij)["abij"];
    //(*H22)["cbjkabij"] -= (*ctfHabij)["acik"];
    //(*H22)["adliabij"] -= (*ctfHabij)["bdjl"];
    allocatedTensorArgument("SimlarityTransformedHamiltonian22", H22);
    allocatedTensorArgument("SimlarityTransformedHamiltonian20", H20);
    allocatedTensorArgument("EnergyDenominators", new CTF::Tensor<>(*ctfDabij));
  }

  determineEnergyShift();

  auto beta(TCC::tensor(std::vector<size_t>(), "beta"));
  auto delta(TCC::tensor(std::vector<size_t>(), "delta"));
  auto energy(TCC::tensor(std::vector<size_t>(), "energy"));
  CTF::Scalar<> ctfBeta(*Cc4s::world);
  CTF::Scalar<> ctfDelta(*Cc4s::world);
  CTF::Scalar<> ctfEnergy(*Cc4s::world);

  DefaultRandomEngine random;
  std::normal_distribution<double> normalDistribution(0.0, 1.0);
  setRandomTensor(*ctfRabij, normalDistribution, random);
  (*ctfLabij)["abij"] = (*ctfRabij)["abij"];

  // Build beta
  ctfBeta[""]  = 0.5*spins*spins*(*ctfLabij)["abij"] * (*ctfRabij)["abij"];
  ctfBeta[""] += -1.0*0.5*spins*(*ctfLabij)["abij"] * (*ctfRabij)["abji"];

  (*ctfLabij)["abij"] = ( 1 / ctfBeta.get_val() ) * (*ctfLabij)["abij"];

  // execute iterations
  int maxIterations(getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS));
  for (int i(0); i < maxIterations; ++i) {
    LOG(0, "DrccdEOM") << "Iteration " << i << "..." << std::endl;

    //buildHbarR
    (*ctfXabij)["abij"] = (*ctfDabij)["abij"] * (*ctfRabij)["abij"];
    (*ctfXabij)["abij"]  += spins * (*ctfRabij)["acik"] * (*ctfHabij)["cbkj"];
    (*ctfXabij)["abij"]  += spins * (*ctfHabij)["caki"] * (*ctfRabij)["bcjk"];
    (*ctfXabij)["abij"]  += -1.0*energyShift * (*ctfRabij)["abij"];

    //buildLHbarR
    ctfEnergy[""] = 0.5*spins*spins * (*ctfXabij)["abij"] * (*ctfLabij)["abij"];
    ctfEnergy[""] += -1.0*0.5 * spins * (*ctfXabij)["abij"] * (*ctfLabij)["abji"];

    //buildNewR
    (*ctfYabij)["abij"] = -1.0*ctfBeta[""] * (*ctfPrevRabij)["abij"];
    (*ctfYabij)["abij"]  += (*ctfXabij)["abij"];
    (*ctfYabij)["abij"]  += -1.0*ctfEnergy[""] * (*ctfRabij)["abij"];
    (*ctfPrevRabij)["abij"] = (*ctfRabij)["abij"];
    (*ctfRabij)["abij"] = (*ctfYabij)["abij"];

    //buildLHbar
    (*ctfXabij)["abij"] = (*ctfVabij)["abij"];
    (*ctfXabij)["abij"]  += (*ctfDabij)["abij"] * (*ctfLabij)["abij"];
    (*ctfXabij)["abij"]  += spins * (*ctfVabij)["acik"] * (*ctfHabij)["bcjk"];
    (*ctfXabij)["abij"]  += spins * (*ctfHabij)["acik"] * (*ctfLabij)["bcjk"];
    (*ctfXabij)["abij"]  += -1.0*energyShift * (*ctfLabij)["abij"];

    //buildNewL
    (*ctfYabij)["abij"] = -1.0*ctfDelta[""] * (*ctfPrevLabij)["abij"];
    (*ctfYabij)["abij"]  += (*ctfXabij)["abij"];
    (*ctfYabij)["abij"]  += -1.0*ctfEnergy[""] * (*ctfLabij)["abij"];
    (*ctfPrevLabij)["abij"] = (*ctfLabij)["abij"];
    (*ctfLabij)["abij"] = (*ctfYabij)["abij"];

    // Build coefficients
    ctfBeta[""] = 0.5*spins*spins*(*ctfLabij)["abij"] * (*ctfRabij)["abij"];
    ctfBeta[""] += -1.0*0.5*spins*(*ctfLabij)["abij"] * (*ctfRabij)["abji"];
    //ctfBeta[""] = beta->getMachineTensor()->tensor[""];
    ctfDelta[""] = std::sqrt(std::abs( ctfBeta.get_val() ));
    ctfBeta[""] = ctfBeta.get_val() / ctfDelta.get_val();


    (*ctfRabij)["abij"] = ( 1 / ctfBeta.get_val() ) * (*ctfRabij)["abij"];
    (*ctfLabij)["abij"] = ( 1 / ctfDelta.get_val() ) * (*ctfLabij)["abij"];

    double b(ctfBeta.get_val()), d(ctfDelta.get_val());
    LOG(1, "DrccdEOM") << "beta=" << b << ", delta=" << d << std::endl;

    double e(ctfEnergy.get_val());

    LOG(1, "DrccdEOM") << "e= " << e+energyShift << std::endl;
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
