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

void DrccdEnergyFromCoulombIntegrals::iterate(int i, Mixer<double> *mixer) {
  iterateTemplate<double>(int i, mixer);
}
void DrccdEnergyFromCoulombIntegrals::iterate(int i, Mixer<complex> *mixer) {
  iterateTemplate<complex>(int i, mixer);
}

template <typename F>
void DrccdEnergyFromCoulombIntegrals::iterateTemplate(int i, Mixer<F> *mixer) {
  Tensor<F> *Tabij(&mixer->getNext());
  // Read Vabij
  Tensor<F> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));

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
    Rabij["abij"] += 2.0 * (*Vabij)["acik"] * (*Tabij)["cbkj"];
    if (linearized) {
      Rabij["abij"] += 2.0 * (*Vabij)["cbkj"] * (*Tabij)["acik"];
    } else {
      // Construct intermediates
      Tensor<F> Cabij(false, *Vabij);
      Cabij["abij"]  = 2.0 * (*Vabij)["cbkj"] * (*Tabij)["acik"];
      Rabij["abij"] += Cabij["abij"];
      Rabij["abij"] += 2.0 * Cabij["acik"] * (*Tabij)["cbkj"];
    }
  }

  // Calculate the amplitdues from the residuum
  doublesAmplitudesFromResiduum(Rabij);
  // And append them to the mixer
  TabijMixer->append(Rabij);
}


void DrccdEnergyFromCoulombIntegrals::dryIterate() {
  {
    // TODO: the Mixer should provide a DryTensor in the future
    // Read the DRCCD amplitudes Tabij
    //DryTensor<> *Tabij(
    getTensorArgument<double, DryTensor<double>>("DrccdDoublesAmplitudes");
    //);

    // Read the Coulomb Integrals Vabij
    DryTensor<> *Vabij(getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals"));
  
    // Allocate Tensors for T2 amplitudes
    DryTensor<> Rabij(*Vabij);

    // Define intermediates
    DryTensor<> Cabij(*Vabij);

    // TODO: implment dryDoublesAmplitudesFromResiduum
    // at the moment, assume usage of Dabij
    DryTensor<> Dabij(*Vabij);
  }
}

