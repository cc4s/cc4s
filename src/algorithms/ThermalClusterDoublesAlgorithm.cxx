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
  recursionLength = getIntegerArgument(
    "recursionLength", DEFAULT_RECURSION_LENGTH
  );
  recursionScale= getRecursionScale(recursionLength);
  LOG(1, getCapitalizedAbbreviation())
    << "Number of imaginary time scales retained=" << recursionLength+1
    << std::endl;
  LOG(1, getCapitalizedAbbreviation())
    << "Imaginary time recursion scaling=" << recursionScale << std::endl;

  energies.resize(recursionLength+1);
  amplitudes.resize(recursionLength+1);
  for (unsigned int i(0); i < energies.size(); ++i) {
    energies[i] = new Scalar<>(*Vabij->wrld);
    amplitudes[i] = new Tensor<>(4, vvoo, syms, *Vabij->wrld, "Tabij");
  }

  // allocate and initialize the eigenenergy difference matrix
  Dai = new Matrix<>(Nv, No, *Vabij->wrld);
  Dai->set_name("Dai");
  Tensor<> *epsi(getTensorArgument<>("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ThermalParticleEigenEnergies"));
  (*Dai)["ai"] =  (*epsa)["a"];
  (*Dai)["ai"] -= (*epsi)["i"];

  beta = 1 / getRealArgument("Temperature");
  LOG(1, getCapitalizedAbbreviation()) << "beta=" << beta << std::endl;

  int n;
  // TODO: determine N from spectrum and temperature
  int N(10);
  double energy;
  for (n = N; n > 0; --n) {
    LOG(0, getCapitalizedAbbreviation()) << "imaginary time scale=" << n << std::endl;
    // update the next amplitudes according to the algorithm implementation
    update(n);
    energy = energies[0]->get_val();
  }

  // TODO: only the sequence of amplitudes at all sampled imaginary times
  // would be of interest.
/*
  std::stringstream amplitudesName;
  amplitudesName << "Thermal" << getAbbreviation() << "DoublesAmplitudes";
  allocatedTensorArgument(
    amplitudesName.str(), Tabij[n&1]
  );
*/
  // currently, we can just dispose them
  for (int i(0); i <= recursionLength; ++i) {
    delete energies[i];
    delete amplitudes[i];
  }
  delete Dai;

  std::stringstream energyName;
  energyName << "Thermal" << getAbbreviation() << "Energy";
  setRealArgument(energyName.str(), energy);
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

  getIntegerArgument("recursionLength", DEFAULT_RECURSION_LENGTH);

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

std::string ThermalClusterDoublesAlgorithm::getCapitalizedAbbreviation() {
  std::string abbreviation(getAbbreviation());
  std::transform(
    abbreviation.begin(), abbreviation.end(), 
    abbreviation.begin(), ::toupper
  );
  return abbreviation;
}

double ThermalClusterDoublesAlgorithm::getRecursionScale(const int M) {
  double q(2.0);
  double delta;
  do {
    const double qM(std::pow(q,M));
    q -= delta = (qM*(q-1) - 1) / (qM*(M*(1-1/q)+1));
  } while (std::abs(delta) > 1e-15);
  return q;
}

