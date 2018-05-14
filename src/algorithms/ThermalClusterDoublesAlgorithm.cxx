#include <algorithms/ThermalClusterDoublesAlgorithm.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <math/IterativePseudoInverse.hpp>
#include <util/LapackMatrix.hpp>
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
  LOG(1, getCapitalizedAbbreviation()) << "beta=" << beta << std::endl;

  // get imaginary time and frequency grids on all nodes
  auto tn( getTensorArgument<real>("ImaginaryTimePoints") );
  std::vector<real> taus(tn->lens[0]);
  tn->read_all(taus.data());
  auto twn( getTensorArgument<real>("ImaginaryTimeWeights") );
  std::vector<real> tauWeights(twn->lens[0]);
  twn->read_all(tauWeights.data());

  auto Vabij( getTensorArgument<real>("PPHHCoulombIntegrals") );
  // get energy differences for propagation
  Dabij = NEW(CTF::Tensor<real>, Vabij->order, Vabij->lens, Vabij->sym);
  fetchDelta(*Dabij);

  // allocate doubles amplitudes on imaginary time grid
  Tabijn.resize(taus.size());
  for (size_t n(0); n < taus.size(); ++n) {
    Tabijn[n] = NEW(CTF::Tensor<real>, false, *Vabij);
  }

  real energy;
  real spins(2.0);
  // number of iterations for determining the amplitudes
  int maxIterationsCount(
    getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS)
  );
  for (int i(0); i < maxIterationsCount; ++i) {
    LOG(0, getCapitalizedAbbreviation()) << "iteration: " << i+1 << std::endl;
    real tau0(0.0);
    CTF::Tensor<real> T0abij(false, *Vabij);
    CTF::Tensor<real> S1abij(false, *Vabij);
    CTF::Scalar<real> direct, exchange;
    for (size_t n(0); n < taus.size(); ++n) {
      real DTau(taus[n] - tau0);
      // energy contribution from previously convolved amplitudes S^I(tau_n-1)
      direct[""] += 0.5*spins*spins* DTau/2 * S1abij["abij"] * (*Vabij)["abij"];
      exchange[""] -= 0.5*spins    * DTau/2 * S1abij["abij"] * (*Vabij)["abji"];

      // get T0abij at tau_n-1 and write previously convolved amplitudes
      // note T(0) is implicitly 0
      if (n > 0) {
        T0abij["abij"] = (*Tabijn[n-1])["abij"];
        (*Tabijn[n-1]) = S1abij["abij"];
      }

      LOG(1, getCapitalizedAbbreviation())
        << "convolving amplitudes at tau_" << n << "=" << taus[n] << std::endl;

      // propagate previously convolved amplitudes S^I(tau_n-1) to this tau_n
      Transform<real, real>(
        std::function<void(real, real &)>( ImaginaryTimePropagation(DTau) )
      ) (
        (*Dabij)["abij"], S1abij["abij"]
      );

      // apply hamiltonian between tau_n-1 and tau_n to update to S^I(tau_n)
      applyHamiltonian(T0abij, *Tabijn[n], DTau, S1abij);

      // energy contribution from convolved amplitudes S^I(tau_n)
      direct[""] += 0.5*spins*spins* DTau/2 * S1abij["abij"] * (*Vabij)["abij"];
      exchange[""] -= 0.5*spins    * DTau/2 * S1abij["abij"] * (*Vabij)["abji"];
    }
    real d(direct.get_val());
    real x(exchange.get_val());
    energy = d + x;
    LOG(2, getCapitalizedAbbreviation()) << "e_d=" << d << std::endl;
    LOG(2, getCapitalizedAbbreviation()) << "e_x=" << x << std::endl;
    LOG(1, getCapitalizedAbbreviation()) << "e=" << energy << std::endl;
  }

  std::stringstream energyName;
  energyName << "Thermal" << getAbbreviation() << "Energy";
  setRealArgument(energyName.str(), energy);
}

void ThermalClusterDoublesAlgorithm::dryRun() {
}

std::string ThermalClusterDoublesAlgorithm::getCapitalizedAbbreviation() {
  std::string abbreviation(getAbbreviation());
  std::transform(
    abbreviation.begin(), abbreviation.end(), 
    abbreviation.begin(), ::toupper
  );
  return abbreviation;
}

std::string ThermalClusterDoublesAlgorithm::getAmplitudeIndices(Tensor<> &T) {
  char indices[T.order+1];
  const int excitationLevel(T.order/2);
  for (int i(0); i < excitationLevel; ++i) {
    indices[i] = static_cast<char>('a'+i);
    indices[i+excitationLevel] = static_cast<char>('i'+i);
  }
  indices[T.order] = 0;
  return indices;
}

void ThermalClusterDoublesAlgorithm::fetchDelta(Tensor<> &Delta) {
  Tensor<> *epsi(getTensorArgument<>("ThermalHoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ThermalParticleEigenEnergies"));
  std::string indices(getAmplitudeIndices(Delta));
  real factor(0.0);
  const int excitationLevel(Delta.order/2);
  for (int i(0); i < excitationLevel; ++i) {
    char aIndex[] = {static_cast<char>('a'+i), 0};
    Delta.sum(+1.0, *epsa,aIndex, factor,indices.c_str());
    factor = 1.0;
    char iIndex[] = {static_cast<char>('i'+i), 0};
    Delta.sum(-1.0, *epsi,iIndex, 1.0,indices.c_str());
  }
}

void ThermalClusterDoublesAlgorithm::thermalContraction(Tensor<> &T) {
  Tensor<> *Ni(getTensorArgument<>("ThermalHoleOccupancies"));
  Tensor<> *Na(getTensorArgument<>("ThermalParticleOccupancies"));
  std::string indices(getAmplitudeIndices(T));
  const int excitationLevel(T.order/2);
  for (int i(0); i < excitationLevel; ++i) {
    char aIndex[] = {static_cast<char>('a'+i), 0};
    T[indices.c_str()] *= (*Na)[aIndex];;
    char iIndex[] = {static_cast<char>('i'+i), 0};
    T[indices.c_str()] *= (*Ni)[iIndex];;
  }
}

