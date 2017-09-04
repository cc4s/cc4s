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
  beta = 1 / getRealArgument("Temperature");

  recursionLength = getIntegerArgument(
    "recursionLength", DEFAULT_RECURSION_LENGTH
  );
  recursionScaling = getRecursionScaling(recursionLength);
  double energyScale( getEnergyScale() );
  int64_t N(std::ceil(std::log(beta*energyScale) / std::log(recursionScaling)));
  N += getIntegerArgument("minIterations", DEFAULT_MIN_ITERATIONS);
  N = std::min(N, getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS));

  LOG(1, getCapitalizedAbbreviation()) << "beta=" << beta << std::endl;
  LOG(1, getCapitalizedAbbreviation())
    << "Imaginart time recursion length=" << recursionLength << std::endl;
  LOG(1, getCapitalizedAbbreviation())
    << "Imaginary time recursion scaling q=" << recursionScaling << std::endl;
  LOG(1, getCapitalizedAbbreviation())
    << "Imaginary time recursion iterations=" << N << std::endl;

  initializeRecursion(N);

  int n;
  double energy;
  for (n = N; n > 0; --n) {
    const double betam( beta*std::pow(recursionScaling,-(n-1)) );
    LOG(0, getCapitalizedAbbreviation()) <<
      "calculating imaginary time scale=beta*q^-" << n-1 <<
      "=" << betam << std::endl;

    // iterate the next amplitudes according to the algorithm implementation
    iterate(n);

    energy = recurse(n);
    LOG(1, getCapitalizedAbbreviation()) << "e=" << energy << std::endl;
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

void ThermalClusterDoublesAlgorithm::initializeRecursion(const int N) {
  // Read the Coulomb Integrals Vabij required for the energy
  Tensor<> *Vabij(getTensorArgument<>("ThermalPPHHCoulombIntegrals"));
  Tensor<> *epsi(getTensorArgument<>("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ThermalParticleEigenEnergies"));
  Tensor<> Dabij(false, *Vabij);

  // Allocate the doubles amplitudes
  energies.resize(recursionLength+1);
  amplitudes.resize(recursionLength+1);
  for (int i(recursionLength); i >= 0; --i) {
    const int n(N+i);
    const double betan( beta*std::pow(recursionScaling,-n) );
    LOG(0, getCapitalizedAbbreviation()) <<
      "initializing imaginary time scale=beta*q^-" << n <<
      "=" << betan << std::endl;
    energies[i] = new Scalar<>(*Vabij->wrld);
    energies[i]->set_name("F");
    amplitudes[i] = new Tensor<>(*Vabij);
    amplitudes[i]->set_name("Tabij");
    Dabij["abij"] =  (*epsa)["a"];
    Dabij["abij"] += (*epsa)["b"];
    Dabij["abij"] -= (*epsi)["i"];
    Dabij["abij"] -= (*epsi)["j"];
    SameSideConnectedImaginaryTimePropagation propagation(betan);
    Dabij.sum(1.0, Dabij,"abij", 0.0,"abij", Univar_Function<>(propagation));
    (*amplitudes[i])["abij"] *= (-1.0) * Dabij["abij"];
  }
}

double ThermalClusterDoublesAlgorithm::recurse(const int n) {
  Tensor<> *Vabij(getTensorArgument<>("ThermalPPHHCoulombIntegrals"));
  Tensor<> *epsi(getTensorArgument<>("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ThermalParticleEigenEnergies"));
  Tensor<> Dabij(false, *Vabij);

  Tensor<> nextT( *amplitudes[recursionLength] );
  Dabij["abij"] =  (*epsa)["a"];
  Dabij["abij"] += (*epsa)["b"];
  Dabij["abij"] -= (*epsi)["i"];
  Dabij["abij"] -= (*epsi)["j"];
  const int m(n+recursionLength);
  const double betam( beta*std::pow(recursionScaling,-m) );
  FreeImaginaryTimePropagation freePropagation(betam);
  Dabij.sum(1.0, Dabij,"abij", 0.0,"abij", Univar_Function<>(freePropagation));
  nextT["abij"] += (*amplitudes[0])["abij"] * Dabij["abij"];

  Dabij["abij"] =  (*epsa)["a"];
  Dabij["abij"] += (*epsa)["b"];
  Dabij["abij"] -= (*epsi)["i"];
  Dabij["abij"] -= (*epsi)["j"];
  SameSideConnectedImaginaryTimePropagation propagation(betam);
  Dabij.sum(1.0, Dabij,"abij", 0.0,"abij", Univar_Function<>(propagation));
  Dabij["abij"] *= (*Vabij)["abij"];

  Tensor<> *Ni(getTensorArgument<>("ThermalHoleOccupancies"));
  Tensor<> *Na(getTensorArgument<>("ThermalParticleOccupancies"));
  Dabij["abij"] *= (*Na)["a"];
  Dabij["abij"] *= (*Na)["b"];
  Dabij["abij"] *= (*Ni)["i"];
  Dabij["abij"] *= (*Ni)["j"];

  Scalar<> nextF( *energies[recursionLength] );
  nextF[""] += (*energies[0])[""];
  nextF[""] -= 2.0 * Dabij["abij"] * (*amplitudes[0])["abij"];
  nextF[""] += 1.0 * Dabij["abji"] * (*amplitudes[0])["abij"];

  // rotate amplitude and energy pointers to make the T[M] the new head T[0]
  Tensor<> *nextTPointer(amplitudes[recursionLength]);
  Scalar<> *nextFPointer(energies[recursionLength]);
  for (int i(recursionLength); i > 0; --i) {
    amplitudes[i] = amplitudes[i-1];
    energies[i] = energies[i-1];
  }
  amplitudes[0] = nextTPointer;
  energies[0] = nextFPointer;
  (*nextTPointer)["abij"] = nextT["abij"];
  (*nextFPointer)[""] = nextF[""];

  return energies[0]->get_val()/beta;
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

