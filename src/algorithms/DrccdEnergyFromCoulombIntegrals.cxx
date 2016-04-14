#include <algorithms/DrccdEnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <util/DryTensor.hpp>
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

void DrccdEnergyFromCoulombIntegrals::iterate(int i) {
  // Read the DRCCD amplitudes Tabij
  Tensor<> *Tabij(&TabijMixer->getNext());
  
  // Read tensors (Vabij, epsi, epsa, and intermediates)
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
  Tensor<> Rabij(false, *Vabij);
  Tensor<> Cabij(false, *Vabij);
  Tensor<> Dabij(false, *Vabij);

  Rabij["abij"] = (*Vabij)["abij"];
  Rabij["abij"] += 2.0 * (*Vabij)["acik"] * (*Tabij)["cbkj"];
  Cabij["abij"] =  2.0 * (*Vabij)["cbkj"] * (*Tabij)["acik"];
  Rabij["abij"] += Cabij["abij"];
  Rabij["abij"] += 2.0 * Cabij["acik"] * (*Tabij)["cbkj"];

  // calculate the amplitdues from the residuum
  doublesAmplitudesFromResiduum(Rabij);
  // and append them to the mixer
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

    // Read the Particle/Hole Eigenenergies epsi epsa
    DryTensor<> *epsi(getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies"));
    DryTensor<> *epsa(getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies"));
  
    // Allocate Tensors for T2 amplitudes
    DryTensor<> Rabij(*Vabij);

    // Define intermediates
    DryTensor<> Cabij(*Vabij);

    // TODO: implment dryDoublesAmplitudesFromResiduum
    // at the moment, assume usage of Dabij
    DryTensor<> Dabij(*Vabij);
  }
}

