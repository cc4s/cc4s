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

  auto Vabij( getTensorArgument<real>("ThermalPPHHCoulombIntegrals") );
  // get energy differences for propagation
  Dabij = NEW(CTF::Tensor<real>, Vabij->order, Vabij->lens, Vabij->sym);
  fetchDelta(*Dabij);

  // number of iterations for determining the amplitudes at each point in time
  int I( getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS) );
  // allocate doubles amplitudes on imaginary time grid
  Tabijn.resize(taus.size());
  for (size_t n(0); n < taus.size(); ++n) {
    Tabijn[n] = NEW(CTF::Tensor<real>, false, *Vabij);
  }
  std::vector<real> energies(taus.size());

  int R( getIntegerArgument("renormalizations", 5) );

  real energy;
  real spins(2.0);
  CTF::Tensor<real> T0abij(false, *Vabij);
  CTF::Tensor<real> S1abij(false, *Vabij);
  real d, x;
  for (int r(-R); r <= 0; ++r) {
    real scale( std::pow(taus.back()/taus.front(),r) );
    LOG(0, getCapitalizedAbbreviation()) << "renormalization level: " << r << std::endl;
    for (size_t N(0); N < taus.size(); ++N) {
      real lastEnergy(0);
      for (int i(0); i < I; ++i) {
        LOG(0, getCapitalizedAbbreviation()) << "iteration: " << i+1 << std::endl;
        real tau0(0.0);
        T0abij["abij"] = 0.0;
        S1abij["abij"] = 0.0;
        CTF::Scalar<real> direct, exchange;
        for (size_t n(0); n <= N; ++n) {
          real tau1(scale*taus[n]);
          real DTau(tau1-tau0);
          // energy contribution from previously convolved amplitudes S^I(tau_n-1)
          direct[""] += 0.5*spins*spins* DTau/2 * S1abij["abij"] * (*Vabij)["abij"];
          exchange[""] -= 0.5*spins    * DTau/2 * S1abij["abij"] * (*Vabij)["abji"];

          // get T0abij at tau_n-1 and write previously convolved amplitudes
          // note T(0) is implicitly 0
          if (n > 0) {
            T0abij["abij"] = (*Tabijn[n-1])["abij"];
            (*Tabijn[n-1])["abij"] = S1abij["abij"];
          }

          LOG(1, getCapitalizedAbbreviation())
            << "convolving amplitudes at tau_" << (n+1) << "=" << tau1
            << std::endl;

          // propagate previously convolved amplitudes S^I(tau_n-1) to this tau_n
          Transform<real, real>(
            std::function<void(real, real &)>( ImaginaryTimePropagation(DTau) )
          ) (
            (*Dabij)["abij"], S1abij["abij"]
          );

          // apply hamiltonian between tau_n-1 and tau_n
          applyHamiltonian(T0abij, *Tabijn[n], DTau, S1abij);

          // energy contribution from convolved amplitudes S^I(tau_n)
          direct[""] += 0.5*spins*spins* DTau/2 * S1abij["abij"] * (*Vabij)["abij"];
          exchange[""] -= 0.5*spins    * DTau/2 * S1abij["abij"] * (*Vabij)["abji"];
          d = direct.get_val();
          x = exchange.get_val();
          real a(S1abij.norm2());
          LOG(2, getCapitalizedAbbreviation()) << "F_d=" << d/tau1 << std::endl;
          LOG(2, getCapitalizedAbbreviation()) << "F_x=" << x/tau1 << std::endl;
          LOG(2, getCapitalizedAbbreviation()) << "|T|=" << a << std::endl;
          tau0 = tau1;
        }
        (*Tabijn[N])["abij"] = S1abij["abij"];
        energies[N] = energy = d + x;
        LOG(1, getCapitalizedAbbreviation()) << "F=" << energy/tau0 << std::endl;
        if (std::abs(1-lastEnergy/energy) < 1e-6) break;
        lastEnergy = energy;
      }
    }
  }

  if (isArgumentGiven("plotAmplitudes")) {
    auto Aai(new CTF::Tensor<real>(2, &S1abij.lens[1], &S1abij.sym[1]));
    (*Aai)["ai"] = S1abij["aaii"];
    allocatedTensorArgument<real>("plotAmplitudes", Aai);
  }
  if (isArgumentGiven("plotDeltas")) {
    auto Dai(new CTF::Tensor<real>(2, &S1abij.lens[1], &S1abij.sym[1]));
    (*Dai)["ai"] = (*Dabij)["aaii"];
    allocatedTensorArgument<real>("plotDeltas", Dai);
  }
  if (isArgumentGiven("plotCoulomb")) {
    auto Vai(new CTF::Tensor<real>(2, &S1abij.lens[1], &S1abij.sym[1]));
    (*Vai)["ai"] = (*Vabij)["aaii"];
    allocatedTensorArgument<real>("plotCoulomb", Vai);
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
    T[indices.c_str()] *= (*Na)[aIndex];
    char iIndex[] = {static_cast<char>('i'+i), 0};
    T[indices.c_str()] *= (*Ni)[iIndex];
  }
}

