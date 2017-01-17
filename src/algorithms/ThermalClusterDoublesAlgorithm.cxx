#include <algorithms/ThermalClusterDoublesAlgorithm.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <Options.hpp>

using namespace CTF;
using namespace cc4s;

ThermalClusterDoublesAlgorithm::ThermalClusterDoublesAlgorithm(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ThermalClusterDoublesAlgorithm::~ThermalClusterDoublesAlgorithm() {
  
}

void ThermalClusterDoublesAlgorithm::run() {
  // Read the Coulomb Integrals Vabij required for the energy
  Tensor<> *Vabij(getTensorArgument<>("ThermalPPHHCoulombIntegrals"));

  // Allocate the doubles amplitudes
  int No(Vabij->lens[2]);
  int Nv(Vabij->lens[0]);
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { Nv, Nv, No, No };
  Tabij[0] = new Tensor<>(4, vvoo, syms, *Vabij->wrld, "Tabij0");
  Tabij[1] = new Tensor<>(4, vvoo, syms, *Vabij->wrld, "Tabij1");

  // Allocate the energy
  directEnergy = new Scalar<>(*Vabij->wrld);
  directEnergy->set_name("directEnergy");
  exchangeEnergy = new Scalar<>(*Vabij->wrld);
  exchangeEnergy->set_name("exchangeEnergy");

  std::string abbreviation(getAbbreviation());
  std::transform(
    abbreviation.begin(), abbreviation.end(), 
    abbreviation.begin(), ::toupper
  );

  // allocate and initialize the eigenenergy difference matrix
  Dai = new Matrix<>(Nv, No, *Vabij->wrld);
  Dai->set_name("Dai");
  Tensor<> *epsi(getTensorArgument<>("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ThermalParticleEigenEnergies"));
  (*Dai)["ai"] =  (*epsa)["a"];
  (*Dai)["ai"] -= (*epsi)["i"];

  beta = 1 / getRealArgument("Temperature");
  LOG(1, abbreviation) << "beta=" << beta << std::endl;

  double dire(0), exce(0);
  samples = getIntegerArgument("samples", DEFAULT_SAMPLES);
  int n;
  for (n = 0; n < samples; ++n) {
    LOG(0, abbreviation) << "step=" << n+1 << std::endl;
    // update the next amplitudes according to the algorithm implementation
    update(n);
    dire = directEnergy->get_val();
    exce = exchangeEnergy->get_val();
    LOG(0, abbreviation) << "e=" << dire+exce << std::endl;
    LOG(1, abbreviation) << "dir=" << dire << std::endl;
    LOG(1, abbreviation) << "exc=" << exce << std::endl;
  }

  std::stringstream amplitudesName;
  // TODO: only the sequence of amplitudes at all sampled imaginary times
  // would be of interest.
/*
  amplitudesName << "Thermal" << getAbbreviation() << "DoublesAmplitudes";
  allocatedTensorArgument(
    amplitudesName.str(), Tabij[n&1]
  );
*/
  // currently, we can just dispose them
  delete Tabij[0]; delete Tabij[1];
  delete Dai;

  std::stringstream energyName;
  energyName << "Thermal" << getAbbreviation() << "Energy";
  setRealArgument(energyName.str(), dire+exce);
}

void ThermalClusterDoublesAlgorithm::dryRun() {
  // Read the Coulomb Integrals Vabij required for the energy
  getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals");

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("ThermalHoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ThermalParticleEigenEnergies")
  );

  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { Nv, Nv, No, No };
  DryTensor<> Tabij0(4, vvoo, syms, SOURCE_LOCATION);
  DryTensor<> Tabij1(4, vvoo, syms, SOURCE_LOCATION);

  // Allocate the energy e
  DryScalar<> energy(SOURCE_LOCATION);

  getIntegerArgument("samples", DEFAULT_SAMPLES);

  // Call the dry iterate of the actual algorithm, which is left open here
  dryUpdate();

  std::stringstream energyName;
  energyName << "Thermal" << getAbbreviation() << "Energy";
  setRealArgument(energyName.str(), 0.0);
}

void ThermalClusterDoublesAlgorithm::dryUpdate() {
  LOG(0, "CluserDoubles") << "Dry run for update not given for Thermal"
    << getAbbreviation() << std::endl;
}

