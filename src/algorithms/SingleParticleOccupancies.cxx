#include <algorithms/SingleParticleOccupancies.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(SingleParticleOccupancies);

SingleParticleOccupancies::SingleParticleOccupancies(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

SingleParticleOccupancies::~SingleParticleOccupancies() {
}

void SingleParticleOccupancies::run() {
  // Read the DRCCD amplitudes Tabij
  Tensor<> *Tabij(getTensorArgument<>("DoublesAmplitudes"));

  int no(Tabij->lens[2]), nv(Tabij->lens[0]);
  // create particle and hole occupancies
  Vector<> *Ni(new Vector<>(no, *Tabij->wrld, "Ni"));
  Vector<> *Na(new Vector<>(nv, *Tabij->wrld, "No"));

  // calculate <T|T> which is <Psi|Psi> in a linearized CC theory
  Scalar<> TT(*Tabij->wrld);
  TT.set_name("TT");
  TT[""] = (*Tabij)["abij"] * (*Tabij)["abij"];

  // calculate <Psi|Na|Psi>
  (*Na)["a"] = +2.0 * (*Tabij)["abij"] * (*Tabij)["abij"];
  // calculate <Psi|Na|Psi> / <Psi|Psi>
  Bivar_Function<> fDivide(&divide<double>);
  Tensor<> Da(false, *Na);
  Da["a"] = TT[""];
  Na->contract(1.0, *Na,"a", Da,"a", 0.0,"a", fDivide);

  // calculate <Psi|Ni|Psi>
  (*Ni)["i"] = -2.0 * (*Tabij)["abij"] * (*Tabij)["abij"];
  // calculate <Psi|Na|Psi> / <Psi|Psi>
  Tensor<> Di(false, *Ni);
  Di["i"] = TT[""];
  Ni->contract(1.0, *Ni,"i", Di,"i", 0.0,"i", fDivide);
  (*Ni)["i"] += 1.0;

  allocatedTensorArgument<>("ParticleOccupancies", Na);
  allocatedTensorArgument<>("HoleOccupancies", Ni);
}

void SingleParticleOccupancies::dryRun() {
  // Read the DRCCD amplitudes Tabij
  DryTensor<> *Tabij(
    getTensorArgument<double, DryTensor<>>("DoublesAmplitudes")
  );

  int no(Tabij->lens[2]), nv(Tabij->lens[0]);
  // create particle and hole occupancies
  DryVector<> *Ni(new DryVector<>(no));
  DryVector<> *Na(new DryVector<>(nv));

  // calculate <T|T> which is <Psi|Psi> in a linearized CC theory
  DryScalar<> TT();

  // calculate <Psi|Na|Psi>
  DryTensor<> Da(*Na);
  DryTensor<> Di(*Ni);

  allocatedTensorArgument<double, DryTensor<>>("ParticleOccupancies", Na);
  allocatedTensorArgument<double, DryTensor<>>("HoleOccupancies", Ni);
}

