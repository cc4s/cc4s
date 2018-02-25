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
  auto nun( getTensorArgument<real>("EvenImaginaryFrequencyPoints") );
  std::vector<real> nus(nun->lens[0]);
  nun->read_all(nus.data());
  auto nuwn( getTensorArgument<real>("EvenImaginaryFrequencyWeights") );
  std::vector<real> nuWeights(nuwn->lens[0]);
  nuwn->read_all(nuWeights.data());

  // get forward and backward Laplace transform on all nodes
  auto ctfCTvn( getTensorArgument<real>("CosineTransform") );
  auto ctfSTvn( getTensorArgument<real>("SineTransform") );
  auto ctfLTvn( new CTF::Tensor<complex>(2, ctfCTvn->lens, ctfCTvn->sym) );
  (*ctfSTvn)["vn"] *= -1;
  toComplexTensor(*ctfCTvn, *ctfSTvn, *ctfLTvn);
  double eta( getRealArgument("EnergyShift", 1.0) );
  // apply shift from forward transformation, which is not done in VASP
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [eta](real tau, complex &T) {
        T *= std::exp(-eta*tau);
      }
    )
  ) (
    (*tn)["n"], (*ctfLTvn)["vn"]
  );
  IterativePseudoInverse<complex> ctfInvLTnv(*ctfLTvn);
  LapackMatrix<complex> LTvn(*ctfLTvn);
  LapackMatrix<complex> invLTnv(ctfInvLTnv.get());

  auto Vabij( getTensorArgument<complex>("PPHHCoulombIntegrals") );

  // test Laplace transform
  CTF::Tensor<real> Dai(2, &Vabij->lens[1], &Vabij->sym[1]);
  fetchDelta(Dai);
  int lens[] = {Vabij->lens[0], Vabij->lens[2], tn->lens[0]};
  int syms[] = {NS, NS, NS};
  auto Tain( new CTF::Tensor<complex>(3, lens, syms) );
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [](real Delta, complex &T) {
        T = Delta;
      }
    )
  ) (
    Dai["ai"], (*Tain)["ain"]
  );
  auto Taiv( new CTF::Tensor<complex>(*Tain) );
  // compute propagator in imaginary time
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [](real tau, complex &T) {
        T = std::exp(-T*tau);
      }
    )
  ) (
    (*tn)["n"], (*Tain)["ain"]
  );
  // compute propagator in imaginary frequency
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [eta](real nu, complex &T) {
        T = 1.0 / complex(std::real(T)+eta, nu);
      }
    )
  ) (
    (*nun)["v"], (*Taiv)["aiv"]
  );
  // numericall transform from time to frequency
  auto NTaiv( new CTF::Tensor<complex>(false, *Taiv) );
  (*NTaiv)["aiv"] = (*Tain)["ain"] * (*ctfLTvn)["vn"];
  CTF::Tensor<complex> diffTaiv(*NTaiv);
  diffTaiv["aiv"] += (-1.0) * (*Taiv)["aiv"];
  double forwardError(frobeniusNorm(diffTaiv));
  LOG(0, getCapitalizedAbbreviation())
    << "forward Laplace transform error=" << forwardError << std::endl;
  // numericall transform from frequency to time
  auto NTain( new CTF::Tensor<complex>(false, *Tain) );
  (*NTain)["aiv"] = (*Taiv)["aiv"] * ctfInvLTnv.get()["vn"];
  CTF::Tensor<complex> diffTain(*NTain);
  diffTain["ain"] += (-1.0) * (*Tain)["ain"];
  double inverseError(frobeniusNorm(diffTain));
  LOG(0, getCapitalizedAbbreviation())
    << "inverse Laplace transform error=" << inverseError << std::endl;

  allocatedTensorArgument<complex>(
    "LaplaceTransform", ctfLTvn
  );
  allocatedTensorArgument<complex>(
    "ExactImaginaryFrequencyPropagators", Taiv
  );
  allocatedTensorArgument<complex>(
    "ExactImaginaryFrequencyPropagators", Taiv
  );
  allocatedTensorArgument<complex>(
    "NumericalImaginaryFrequencyPropagators", NTaiv
  );

  allocatedTensorArgument<complex>(
    "ExactImaginaryTimePropagators", Tain
  );
  allocatedTensorArgument<complex>(
    "NumericalImaginaryTimePropagators", NTain
  );

  // get energy differences for propagation
  Dabij = NEW(CTF::Tensor<real>, Vabij->order, Vabij->lens, Vabij->sym);
  fetchDelta(*Dabij);

  // allocate doubles amplitudes on imaginary time and frequency grid
  Tabijn.resize(taus.size());
  for (size_t n(0); n < taus.size(); ++n) {
    Tabijn[n] = NEW(CTF::Tensor<complex>, false, *Vabij);
  }
  Tabijv.resize(nus.size());
  for (size_t v(0); v < nus.size(); ++v) {
    Tabijv[v] = NEW(CTF::Tensor<complex>, false, *Vabij);
  }

  double energy;
  // number of iterations for determining the amplitudes
  int maxIterationsCount(
    getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS)
  );
  for (int i(0); i < maxIterationsCount; ++i) {
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

    energy = 0.0;
    LOG(1, getCapitalizedAbbreviation()) << "e=" << energy << std::endl;
  }

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
  double factor(0.0);
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

