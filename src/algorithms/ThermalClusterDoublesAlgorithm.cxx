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
  recursionLength = getIntegerArgument(
    "recursionLength", DEFAULT_RECURSION_LENGTH
  );
  recursionScaling = getRecursionScaling(recursionLength);
  beta = 1 / getRealArgument("Temperature");

  double energyScale( getEnergyScale() );
  int N( std::ceil(std::log(beta*energyScale) / std::log(recursionScaling)) );
  N += getIntegerArgument("minIterations", DEFAULT_MIN_ITERATIONS);

  LOG(1, getCapitalizedAbbreviation()) << "beta=" << beta << std::endl;
  LOG(1, getCapitalizedAbbreviation())
    << "Imaginart time recursion length=" << recursionLength << std::endl;
  LOG(1, getCapitalizedAbbreviation())
    << "Imaginary time recursion scaling=" << recursionScaling << std::endl;
  LOG(1, getCapitalizedAbbreviation())
    << "Imaginary time recursion iterations=" << N << std::endl;


  initializeRecursion();


  int n;
  double energy;
  for (n = N; n > 0; --n) {
    LOG(0, getCapitalizedAbbreviation()) << "imaginary time scale=" << n << std::endl;
    // iterate the next amplitudes according to the algorithm implementation
    
    iterate(n);
    energy = energies[0]->get_val();
  }

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
  dryIterate();

  std::stringstream energyName;
  energyName << "Thermal" << getAbbreviation() << "Energy";
  setRealArgument(energyName.str(), 0.0);
}

void ThermalClusterDoublesAlgorithm::dryIterate() {
  LOG(0, "CluserDoubles") << "Dry run for iterate not given for Thermal"
    << getAbbreviation() << std::endl;
}

void ThermalClusterDoublesAlgorithm::initializeRecursion() {
  // Read the Coulomb Integrals Vabij required for the energy
  Tensor<> *Vabij(getTensorArgument<>("ThermalPPHHCoulombIntegrals"));

  // Allocate the doubles amplitudes
  int No(Vabij->lens[2]);
  int Nv(Vabij->lens[0]);
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { Nv, Nv, No, No };
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
	}

std::string ThermalClusterDoublesAlgorithm::getCapitalizedAbbreviation() {
  std::string abbreviation(getAbbreviation());
  std::transform(
    abbreviation.begin(), abbreviation.end(), 
    abbreviation.begin(), ::toupper
  );
  return abbreviation;
}

/**
 * \brief Calculates and returns the scaling factor \f$q\f$ satisfying
 * \f$q^{M+1}-q^M-1 = 0\f$.
 **/
double ThermalClusterDoublesAlgorithm::getRecursionScaling(const int M) {
  double q(2.0);
  double delta;
  do {
    const double qM(std::pow(q,M));
    q -= delta = (qM*(q-1) - 1) / (qM*(M*(q-1)/q+1));
  } while (std::abs(delta) > 1e-15);
  return q;
}

double ThermalClusterDoublesAlgorithm::getEnergyScale() {
  Tensor<> *epsi(getTensorArgument<>("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ThermalParticleEigenEnergies"));
  std::vector<double> holeEnergies(epsi->lens[0]);
  epsi->read_all(holeEnergies.data());
  std::vector<double> particleEnergies(epsa->lens[0]);
  epsa->read_all(particleEnergies.data());
  return particleEnergies.back() - holeEnergies.front();
}

