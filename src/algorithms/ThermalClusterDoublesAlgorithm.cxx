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

  real eta( getRealArgument("EnergyShift", 1.0) );

  // forward and inverse Laplace transform weights
  auto ctfCTvn( getTensorArgument<real>("CosineTransform") );
  auto ctfSTvn( getTensorArgument<real>("SineTransform") );
  auto ctfICTnv( getTensorArgument<real>("InverseCosineTransform") );
  auto ctfISTnv( getTensorArgument<real>("InverseSineTransform") );

  // for testing
  CTF::Tensor<complex> cCTvn(2, ctfCTvn->lens, ctfCTvn->sym);
  toComplexTensor(*ctfCTvn, cCTvn);
  CTF::Tensor<complex> cSTvn(2, ctfSTvn->lens, ctfSTvn->sym);
  toComplexTensor(*ctfSTvn, cSTvn);
  CTF::Tensor<complex> cICTnv(2, ctfICTnv->lens, ctfICTnv->sym);
  toComplexTensor(*ctfICTnv, cICTnv);
  CTF::Tensor<complex> cISTnv(2, ctfISTnv->lens, ctfISTnv->sym);
  toComplexTensor(*ctfISTnv, cISTnv);

  LapackMatrix<real> CTvn(*ctfCTvn);
  LapackMatrix<real> STvn(*ctfSTvn);
  LapackMatrix<real> ICTnv(*ctfICTnv);
  LapackMatrix<real> ISTnv(*ctfISTnv);

  // cos(nu_v*tau_n) and sin(nu_v*tau_n)
  LapackMatrix<real> Cvn(CTvn);
  LapackMatrix<real> Svn(CTvn);
  for (size_t n(0); n < taus.size(); ++n) {
    for (size_t v(0); v < nus.size(); ++v) {
      Cvn(v,n) = std::cos(nus[v]*taus[n]);
      Svn(v,n) = std::sin(nus[v]*taus[n]);
    }
  }


  // test Laplace transform
  auto ctfCvn( *ctfCTvn );
  Cvn.write(ctfCvn);
  auto ctfSvn( *ctfSTvn );
  Svn.write(ctfSvn);
  CTF::Tensor<complex> cCvn(2, ctfCvn.lens, ctfCvn.sym);
  toComplexTensor(ctfCvn, cCvn);
  CTF::Tensor<complex> cSvn(2, ctfSvn.lens, ctfSvn.sym);
  toComplexTensor(ctfSvn, cSvn);

  auto Vabij( getTensorArgument<complex>("PPHHCoulombIntegrals") );
  CTF::Tensor<real> Dai(2, &Vabij->lens[1], &Vabij->sym[1]);
  fetchDelta(Dai);
  int lens[] = {Vabij->lens[0], Vabij->lens[2], tn->lens[0]};
  int syms[] = {NS, NS, NS};
  auto Tain( new CTF::Tensor<complex>(3, lens, syms) );
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [](real Delta, complex &T) { T = Delta; }
    )
  ) (
    Dai["ai"], (*Tain)["ain"]
  );
  auto TRaiv( new CTF::Tensor<complex>(*Tain) );
  auto TIaiv( new CTF::Tensor<complex>(*Tain) );
  // propagator in imaginary time
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [](real tau, complex &T) { T = std::exp(-std::real(T)*tau); }
    )
  ) (
    (*tn)["n"], (*Tain)["ain"]
  );
  // propagator in imaginary frequency, real part
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [eta](real nu, complex &T) {
        real Delta_( std::real(T)+eta );
        T = Delta_ / (Delta_*Delta_ + nu*nu);
      }
    )
  ) (
    (*nun)["v"], (*TRaiv)["aiv"]
  );
  // propagator in imaginary frequency, imaginary part
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [eta](real nu, complex &T) {
        real Delta_( std::real(T)+eta );
        T = -nu / (Delta_*Delta_ + nu*nu);
      }
    )
  ) (
    (*nun)["v"], (*TIaiv)["aiv"]
  );
  // numericall transform from time to frequency
  auto T_ain( *Tain );
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [eta](real tau, complex &T) { T *= std::exp(-eta*tau); }
    )
  ) (
    (*tn)["n"], T_ain["ain"]
  );
  auto NTRaiv( new CTF::Tensor<complex>(false, *TRaiv) );
  (*NTRaiv)["aiv"] = (+1.0) * T_ain["ain"] * cCTvn["vn"] * cCvn["vn"];
  auto NTIaiv( new CTF::Tensor<complex>(false, *TIaiv) );
  (*NTIaiv)["aiv"] = (-1.0) * T_ain["ain"] * cSTvn["vn"] * cSvn["vn"];
  {
    CTF::Tensor<complex> diffTaiv(*NTRaiv);
    diffTaiv["aiv"] += (-1.0) * (*TRaiv)["aiv"];
    real forwardError(frobeniusNorm(diffTaiv));
    LOG(0, getCapitalizedAbbreviation())
      << "forward Laplace transform error real=" << forwardError << std::endl;
  }
  {
    CTF::Tensor<complex> diffTaiv(*NTIaiv);
    diffTaiv["aiv"] += (-1.0) * (*TIaiv)["aiv"];
    real forwardError(frobeniusNorm(diffTaiv));
    LOG(0, getCapitalizedAbbreviation())
      << "forward Laplace transform error imag=" << forwardError << std::endl;
  }
  allocatedTensorArgument<complex>(
    "ExactImaginaryFrequencyPropagatorsReal", TRaiv
  );
  allocatedTensorArgument<complex>(
    "ExactImaginaryFrequencyPropagatorsImag", TIaiv
  );
  allocatedTensorArgument<complex>(
    "NumericalImaginaryFrequencyPropagatorsReal", NTRaiv
  );
  allocatedTensorArgument<complex>(
    "NumericalImaginaryFrequencyPropagatorsImag", NTIaiv
  );

  // numericall transform from frequency to time
  auto NRTain( new CTF::Tensor<complex>(false, *Tain) );
  (*NRTain)["ain"] = (+4.0) * (*TRaiv)["aiv"] * cICTnv["nv"] * cCvn["vn"];
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [eta](real tau, complex &T) { T *= std::exp(+eta*tau); }
    )
  ) (
    (*tn)["n"], (*NRTain)["ain"]
  );
  auto NITain( new CTF::Tensor<complex>(false, *Tain) );
  (*NITain)["ain"] = (-4.0) * (*TIaiv)["aiv"] * cISTnv["nv"] * cSvn["vn"];
  CTF::Transform<real, complex>(
    std::function<void(real, complex &)>(
      [eta](real tau, complex &T) { T *= std::exp(+eta*tau); }
    )
  ) (
    (*tn)["n"], (*NITain)["ain"]
  );
  {
    CTF::Tensor<complex> diffTain(*NRTain);
    diffTain["ain"] += (-1.0) * (*Tain)["ain"];
    real inverseError(frobeniusNorm(diffTain));
    LOG(0, getCapitalizedAbbreviation())
      << "inverse Laplace transform error real=" << inverseError << std::endl;
  }
  {
    CTF::Tensor<complex> diffTain(*NITain);
    diffTain["ain"] += (-1.0) * (*Tain)["ain"];
    real inverseError(frobeniusNorm(diffTain));
    LOG(0, getCapitalizedAbbreviation())
      << "inverse Laplace transform error imag=" << inverseError << std::endl;
  }
  allocatedTensorArgument<complex>(
    "ExactImaginaryTimePropagators", Tain
  );
  allocatedTensorArgument<complex>(
    "NumericalImaginaryTimePropagatorsReal", NRTain
  );
  allocatedTensorArgument<complex>(
    "NumericalImaginaryTimePropagatorsImag", NITain
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
  // Read the Coulomb Integrals Vabij required for the energy
  getTensorArgument<real, DryTensor<real>>("PPHHCoulombIntegrals");

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(
    getTensorArgument<real, DryTensor<real>>("ThermalHoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<real, DryTensor<real>>("ThermalParticleEigenEnergies")
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

