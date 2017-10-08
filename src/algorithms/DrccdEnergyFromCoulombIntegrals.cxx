#include <algorithms/DrccdEnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DrccdEnergyFromCoulombIntegrals);

DrccdEnergyFromCoulombIntegrals::DrccdEnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): ClusterSinglesDoublesAlgorithm(argumentList) {
}

DrccdEnergyFromCoulombIntegrals::~DrccdEnergyFromCoulombIntegrals() {
}

PTR(FockVector<double>) DrccdEnergyFromCoulombIntegrals::getResiduum(
  const int iteration, const PTR(FockVector<double>) &amplitudes
) {
  return getResiduum<double>(iteration, amplitudes);
}

PTR(FockVector<complex>) DrccdEnergyFromCoulombIntegrals::getResiduum(
  const int iteration, const PTR(FockVector<complex>) &amplitudes
) {
  return getResiduum<complex>(iteration, amplitudes);
}

template <typename F>
PTR(FockVector<F>) DrccdEnergyFromCoulombIntegrals::getResiduum(
  const int iteration, const PTR(FockVector<F>) &amplitudes
) {
  // read all required integrals
  auto Vabij(getTensorArgument<F>("PPHHCoulombIntegrals"));
  auto Vaijb(getTensorArgument<F>("PHHPCoulombIntegrals"));
  auto Vijab(getTensorArgument<F>("HHPPCoulombIntegrals"));

  // get amplitude parts
  auto Tabij( amplitudes->get(1) );

  // construct residuum
  auto residuum( NEW(FockVector<F>, *amplitudes) );
  *residuum *= F(0);
  auto Rabij( residuum->get(1) );

  int linearized(getIntegerArgument("linearized", 0));
  if (linearized) {
    LOG(1, getCapitalizedAbbreviation()) <<
      "Solving linearized T2 Amplitude Equations" << std::endl;
  } else {
    LOG(1, getCapitalizedAbbreviation()) <<
      "Solving T2 Amplitude Equations" << std::endl;
  }

  (*Rabij)["abij"] += (*Vabij)["abij"];

  if (iteration > 0 || isArgumentGiven("startingDoublesAmplitudes")) {
    // for the rest iterations compute the drCCD residuum
    (*Rabij)["abij"] += (*Vabij)["abij"];
    (*Rabij)["abij"] += 2.0 * (*Vaijb)["akic"] * (*Tabij)["cbkj"];
    (*Rabij)["abij"] += 2.0 * (*Vaijb)["bkjc"] * (*Tabij)["acik"];
    if (!linearized) {
      // Construct intermediates
      Tensor<F> Calid(false, *Vaijb);
      Calid["alid"]  = 2.0 * (*Vijab)["klcd"] * (*Tabij)["acik"];
      (*Rabij)["abij"] += 2.0 * Calid["alid"] * (*Tabij)["dblj"];
    }
  }

  return residuum;
}

