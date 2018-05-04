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

  auto Vabij( getTensorArgument<complex>("PPHHCoulombIntegrals") );
  // get energy differences for propagation
  Dabij = NEW(CTF::Tensor<real>, Vabij->order, Vabij->lens, Vabij->sym);
  fetchDelta(*Dabij);

  // allocate doubles amplitudes on imaginary time and frequency grid
  Tabijn.resize(taus.size());
  for (size_t n(0); n < taus.size(); ++n) {
    Tabijn[n] = NEW(CTF::Tensor<real>, false, *Vabij);
  }

  real energy;
  // number of iterations for determining the amplitudes
  int maxIterationsCount(
    getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS)
  );
  for (int i(0); i < maxIterationsCount; ++i) {
/*
    LOG(0, getCapitalizedAbbreviation()) << "iteration: " << i+1 << std::endl;
    for (size_t n(0); n < taus.size(); ++n) {
      LOG(1, getCapitalizedAbbreviation())
        << "computing residues at tau_" << n << "=" << taus[n] << std::endl;
      // iterate amplitudes on time grid
      getResiduum(*Tabijn[n]);
    }
    // transform to frequency grid
    LOG(1, getCapitalizedAbbreviation())
      << "transforming to frequency grid" << std::endl;
    for (size_t v(0); v < nus.size(); ++v) {
      Tabijv[v]->sum(LTvn(v,0), *Tabijn[0],"abij", 0.0,"abij");
      for (size_t n(1); n < taus.size(); ++n) {
        Tabijv[v]->sum(LTvn(v,n), *Tabijn[n],"abij", 1.0,"abij");
      }
    }
    // convolve with propagator
    class Convolution {
    public:
      Convolution(real nu_, real eta_): nu(nu_), eta(eta_) { }
      void operator()(real Delta, complex &T) {
        T /= complex(Delta+eta, nu);
      }
    protected:
      real nu, eta;
    };
    for (size_t v(0); v < nus.size(); ++v) {
      LOG(1, getCapitalizedAbbreviation())
        << "computing convolution at nu_" << v << "=" << nus[v] << std::endl;
      Convolution convolution(nus[v], eta);
      CTF::Transform<real, complex>(
        std::function<void(real, complex &)>(convolution)
      ) (
        (*Dabij)["abij"], (*Tabijv[v])["abij"]
      );
    }

    // transform to time grid
    LOG(1, getCapitalizedAbbreviation())
      << "transforming back to time grid" << std::endl;
    for (size_t n(1); n < taus.size(); ++n) {
      Tabijn[n]->sum(invLTnv(n,0), *Tabijv[0],"abij", 0.0,"abij");
      for (size_t v(0); v < nus.size(); ++v) {
        Tabijn[n]->sum(invLTnv(n,v), *Tabijv[v],"abij", 1.0,"abij");
      }
    }
*/
    energy = 0.0;
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

