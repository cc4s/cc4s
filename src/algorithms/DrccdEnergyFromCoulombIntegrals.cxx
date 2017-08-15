#include <algorithms/DrccdEnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DrccdEnergyFromCoulombIntegrals);

DrccdEnergyFromCoulombIntegrals::DrccdEnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): ClusterDoublesAlgorithm(argumentList) {
}

DrccdEnergyFromCoulombIntegrals::~DrccdEnergyFromCoulombIntegrals() {
}

void DrccdEnergyFromCoulombIntegrals::iterate(
  int i, Mixer<double> *TaiMixer, Mixer<double> *TabijMixer
) {
  iterate<double>(i, TabijMixer);
}
void DrccdEnergyFromCoulombIntegrals::iterate(
  int i, Mixer<complex> *TaiMixer, Mixer<complex> *TabijMixer
) {
  iterate<complex>(i, TabijMixer);
}

template <typename F>
void DrccdEnergyFromCoulombIntegrals::iterate(
  int i, Mixer<F> *TabijMixer
) {
  Tensor<F> *Tabij(&TabijMixer->getNext());
  // Read all required integrals
  Tensor<F> *Vabij(getTensorArgument<F>("PPHHCoulombIntegrals"));
  Tensor<F> *Vaijb(getTensorArgument<F>("PHHPCoulombIntegrals"));
  Tensor<F> *Vijab(getTensorArgument<F>("HHPPCoulombIntegrals"));

  // Construct intermediate Amplitudes
  Tensor<F> Rabij(false, *Tabij);

  std::string abbreviation(getAbbreviation());
  std::transform(abbreviation.begin(), abbreviation.end(), 
                 abbreviation.begin(), ::toupper);

  int linearized(getIntegerArgument("linearized", 0));
  if (linearized) {
    LOG(1, abbreviation) << "Solving linearized T2 Amplitude Equations" << std::endl;
  } else {
    LOG(1, abbreviation) << "Solving T2 Amplitude Equations" << std::endl;
  }

  if (i == 0) {
    // For first iteration compute only the MP2 amplitudes 
    // Since Tabij = 0, Vabij is the only non-zero term

    Rabij["abij"] += (*Vabij)["abij"];
  } 
  else {
    // For the rest iterations compute the DRCCD amplitudes
    Rabij["abij"]  = (*Vabij)["abij"];
    Rabij["abij"] += 2.0 * (*Vaijb)["akic"] * (*Tabij)["cbkj"];
    Rabij["abij"] += 2.0 * (*Vaijb)["bkjc"] * (*Tabij)["acik"];
    if (!linearized) {
      // Construct intermediates
      Tensor<F> Calid(false, *Vaijb);
      Calid["alid"]  = 2.0 * (*Vijab)["klcd"] * (*Tabij)["acik"];
      Rabij["abij"] += 2.0 * Calid["alid"] * (*Tabij)["dblj"];
    }
  }

  // Calculate the amplitdues from the residuum
  doublesAmplitudesFromResiduum(Rabij);
  // And append them to the mixer
  TabijMixer->append(Rabij);
}

