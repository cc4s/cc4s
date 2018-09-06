#include <algorithms/ThermalClusterDoublesAlgorithm.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <math/IterativePseudoInverse.hpp>
#include <util/LapackMatrix.hpp>
#include <util/BlacsWorld.hpp>
#include <util/ScaLapackMatrix.hpp>
#include <util/ScaLapackHermitianEigenSystemDc.hpp>
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

  computeSqrtOccupancies();
  diagonalizeSinglesHamiltonian();

  real tda( getTammDancoffEnergy() );
  LOG(1, getCapitalizedAbbreviation()) << "TDA F=" << tda << std::endl;
  if (isArgumentGiven("ThermalTdaEnergy")) {
    setRealArgument("ThermalTdaEnergy", tda);
  }

  // compute the other contributions perturbatively
  // get imaginary time and frequency grids on all nodes
  auto tn( getTensorArgument<real>("ImaginaryTimePoints") );
  std::vector<real> taus(tn->lens[0]);
  tn->read_all(taus.data());
  // scale to [0,beta]
  for (auto &tau: taus) { tau *= beta; }

  auto Vabij( getTensorArgument<real>("ThermalPPHHCoulombIntegrals") );

  // number of iterations for determining the amplitudes at each point in time
  int I( getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS) );
  // allocate doubles amplitudes on imaginary time grid
  // the source amplitudes are stored in orbital space
  // with fully closed contraction weights, suitable for interpolation
  Tabijn.resize(taus.size());
  for (size_t n(0); n < taus.size(); ++n) {
    Tabijn[n] = NEW(CTF::Tensor<real>, false, *Vabij);
  }
  // the target amplitudes are stored in singles-mode space
  // with half-open contraction weights, suitable for propagation
  std::vector<PTR(Tensor<real>)> SFn(taus.size()), SFGn(taus.size());
  for (size_t n(0); n < taus.size(); ++n) {
    SFGn[n] = NEW(CTF::Tensor<real>, false, *VdFG);
  }
  std::vector<real> energies(taus.size());

  // TODO: test then accept/dismiss renormalizations, for now disabled
  int R( getIntegerArgument("renormalizations", 0) );

  real energy;
  real accuracy(getRealArgument("accuracy", 1e-7));
  CTF::Tensor<real> T0abij(false, *Vabij);
  CTF::Tensor<real> S1FG(false, *VdFG);
  for (int r(-R); r <= 0; ++r) {
    real scale( std::pow(taus.back()/taus.front(),r) );
    LOG(0, getCapitalizedAbbreviation()) << "renormalization level: " << r << std::endl;
    {
      size_t N(taus.size()-1);
//    for (size_t N(0); N < taus.size(); ++N) {
      real lastEnergy(0);
      for (int i(0); i < I; ++i) {
        LOG(0, getCapitalizedAbbreviation()) << "iteration: " << i+1 << std::endl;
        real tau0(0.0);
        T0abij["abij"] = 0.0;
        S1FG["FG"] = 0.0;
        real direct(0.0), exchange(0.0);
        for (size_t n(0); n <= N; ++n) {
          real tau1(scale*taus[n]);
          real DTau(tau1-tau0);
          // energy contribution from previously convolved amplitudes ST^I(tau_n-1)
          computeEnergyContribution(S1FG, DTau, direct, exchange);

          // get T0abij at tau_n-1 and write previously convolved amplitudes
          // note T(0) is implicitly 0
          if (n > 0) {
            T0abij["abij"] = (*Tabijn[n-1])["abij"];
            (*SFGn[n-1])["FG"] = S1FG["FG"];
          }

          LOG(1, getCapitalizedAbbreviation())
            << "convolving amplitudes at tau_" << (n+1) << "=" << tau1
            << std::endl;

          // propagate previously convolved amplitudes ST^I(tau_n-1) to this tau_n
          ImaginaryTimePropagation imaginaryTimePropagation(DTau);
          Transform<real, real>(
            std::function<void(real, real &)>( imaginaryTimePropagation )
          ) (
            (*lambdaFG)["FG"], S1FG["FG"]
          );

          // apply hamiltonian between tau_n-1 and tau_n
          applyHamiltonian(T0abij, *Tabijn[n], DTau, S1FG);

          // energy contribution from convolved amplitudes ST^I(tau_n)
          computeEnergyContribution(S1FG, DTau, direct, exchange);
          LOG(2, getCapitalizedAbbreviation())
            << "F_d=" << direct/tau1 << std::endl;
          LOG(2, getCapitalizedAbbreviation())
            << "F_x=" << exchange/tau1 << std::endl;
          tau0 = tau1;
        }
        (*SFGn[N])["FG"] = S1FG["FG"];
        energies[N] = energy = direct + exchange;
        LOG(1, getCapitalizedAbbreviation()) << "F=" << energy/tau0 << std::endl;
        if (std::abs(1-lastEnergy/energy) < accuracy) break;
        lastEnergy = energy;
        real mixingRatio( getRealArgument("mixingRatio", 1.0) );
        for (size_t n(0); n <= N; ++n) {
          (*Tabijn[n])["abij"] *= (1-mixingRatio);
          // convert to orbital base, fully close contraction weights and mix
          (*Tabijn[n])["abij"] += mixingRatio *
            (*ga)["a"] * (*ga)["b"] * (*gi)["i"] * (*gi)["j"] *
            (*UaiF)["aiF"] * (*UaiF)["bjG"] * (*SFGn[n])["FG"];
          // write norm
          real a(Tabijn[n]->norm2());
          LOG(2, getCapitalizedAbbreviation()) << "|T|=" << a << std::endl;
        }
      }
    }
  }

  std::stringstream energyName;
  energyName << "Thermal" << getAbbreviation() << "Energy";
  setRealArgument(energyName.str(), energy/beta);

  if (isArgumentGiven("plotAmplitudes")) {
    auto newTFG(new CTF::Tensor<real>(S1FG));
    allocatedTensorArgument<real>("plotAmplitudes", newTFG);
  }
  if (isArgumentGiven("plotLambdas")) {
    auto newLambdaFG(new CTF::Tensor<real>(*lambdaFG));
    allocatedTensorArgument<real>("plotLambdas", newLambdaFG);
  }
  if (isArgumentGiven("plotCoulomb")) {
    auto newVdFG(new CTF::Tensor<real>(*VdFG));
    allocatedTensorArgument<real>("plotCoulomb", newVdFG);
  }
  // diagonalizeDoublesAmplitudes();
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

cc4s::real ThermalClusterDoublesAlgorithm::getTammDancoffEnergy() {
  auto Vabij( getTensorArgument<>("ThermalPPHHCoulombIntegrals") );
  real spins( getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0 );

  // compute Tamm-Dancoff Approximation (TDA)
  Scalar<> e;
  int NF(lambdaF->lens[0]);
  // doubles: two particle/hole pairs F&G
  VdFG = NEW(Tensor<real>, 2, std::vector<int>({NF,NF}).data());
  (*VdFG)["FG"] = (*UaiF)["ckF"] * (*UaiF)["dlG"]
    * (*ga)["c"] * (*ga)["d"] * (*gi)["k"] * (*gi)["l"] * (*Vabij)["cdkl"];
  VxFG = NEW(Tensor<real>, false, *VdFG);
  (*VxFG)["FG"] = (*UaiF)["ckF"] * (*UaiF)["dlG"]
    * (*ga)["c"] * (*ga)["d"] * (*gi)["k"] * (*gi)["l"] * (*Vabij)["cdlk"];

  // propagate doubles
  lambdaFG = NEW(Tensor<real>, false, *VdFG);
  (*lambdaFG)["FG"] =  (*lambdaF)["F"];
  (*lambdaFG)["FG"] += (*lambdaF)["G"];
  Tensor<real> TFG(*VdFG);
  SecondOrderIntegral secondOrderIntegral(beta);
  Transform<real, real>(std::function<void(real,real&)>(secondOrderIntegral))(
    (*lambdaFG)["FG"], TFG["FG"]
  );

  e[""] = -0.5 * spins*spins * TFG["FG"] * (*VdFG)["FG"];
  real tdad(e.get_val());
  LOG(1, getCapitalizedAbbreviation()) << "TDA F_d=" << tdad << std::endl;

  e[""] = +0.5 * spins * TFG["FG"] * (*VxFG)["FG"];
  real tdax(e.get_val());
  LOG(1, getCapitalizedAbbreviation()) << "TDA F_x=" << tdax << std::endl;

  if (isArgumentGiven("SinglesHamiltonianWeights")) {
    auto singlesWeights(new Tensor<real>(1, std::vector<int>({NF}).data()) );
    (*singlesWeights)["F"] = 0.5 * spins*spins * TFG["FF"] * (*VdFG)["FF"];
    allocatedTensorArgument<real>("SinglesHamiltonianWeights", singlesWeights);
  }

  return tdad+tdax;
}

void ThermalClusterDoublesAlgorithm::computeEnergyContribution(
  Tensor<real> &SFG, const real DTau,
  real &direct, real &exchange
) {
  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Scalar<real> energy;

  // direct term
  energy[""] = 0.5 * spins * spins * DTau/2 * SFG["FG"] * (*VdFG)["FG"];
  direct += energy.get_val();
  // exchange term
  energy[""] = (-0.5) * spins * DTau/2 * SFG["FG"] * (*VxFG)["FG"];
  exchange += energy.get_val();
}

void ThermalClusterDoublesAlgorithm::computeSqrtOccupancies() {
  // compute g^a and g_i = sqrt(occupancies)
  auto Ni( getTensorArgument<>("ThermalHoleOccupancies") );
  auto Na( getTensorArgument<>("ThermalParticleOccupancies") );
  gi = NEW(Tensor<real>, *Ni);
  Transform<real>(
    std::function<void(real &)>([](real &f) { f = std::sqrt(f); } )
  ) (
    (*gi)["i"]
  );
  ga = NEW(Tensor<real>, *Na);
  Transform<real>(
    std::function<void(real &)>([](real &f) { f = std::sqrt(f); } )
  ) (
    (*ga)["a"]
  );
}

void ThermalClusterDoublesAlgorithm::diagonalizeSinglesHamiltonian() {
  auto epsi(getTensorArgument<>("ThermalHoleEigenEnergies"));
  auto epsa(getTensorArgument<>("ThermalParticleEigenEnergies"));
  auto Vbija(getTensorArgument<>("ThermalPHHPCoulombIntegrals"));
  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);

  // build to Hbjai
  int Nv(epsa->lens[0]); int No(epsi->lens[0]);
  int lens[] = { Nv,No, Nv,No };
  auto Hbjai(NEW(Tensor<>, 4, lens, Vbija->sym, *Vbija->wrld, "Hbjai"));
  // bubble from H_1 has contraction weights
  (*Hbjai)["bjai"] = spins * (*Vbija)["bija"];
  (*Hbjai)["bjai"] *= (*ga)["a"];
  (*Hbjai)["bjai"] *= (*ga)["b"];
  (*Hbjai)["bjai"] *= (*gi)["i"];
  (*Hbjai)["bjai"] *= (*gi)["j"];
  // unperturbed propatation is diagonal
  (*Hbjai)["bjbj"] += (*epsa)["b"];
  (*Hbjai)["bjbj"] -= (*epsi)["j"];

  LOG(1, getCapitalizedAbbreviation())
    << "diagonalizing singles part of Hamiltonian..." << std::endl;
  // H(bj)(ai) = U.S.U^T, seen as a matrix with compound indices
  BlacsWorld world(Hbjai->wrld->rank, Hbjai->wrld->np);
  int NvNo(Nv*No);
  int scaHLens[2] = { NvNo, NvNo };
  auto scaHbjai(NEW(ScaLapackMatrix<>, *Hbjai, scaHLens, &world));
  // release unneeded resources early
  Hbjai = nullptr;

  // use ScaLapack routines to diagonalise the matrix U.Lambda.U^T
  auto scaU(NEW(ScaLapackMatrix<>, *scaHbjai));
  ScaLapackHermitianEigenSystemDc<> eigenSystem(scaHbjai, scaU);
  std::vector<real> lambdas(NvNo);
  eigenSystem.solve(lambdas.data());
  scaHbjai = nullptr;

  // write unitary matrix U(ai)(F) back to CTF as tensor UaiF
  int ULens[3] = { Nv, No, NvNo };
  UaiF = NEW(Tensor<>, 3, ULens, Vbija->sym, *Vbija->wrld, "UaiF");
  scaU->write(*UaiF);
  scaU = nullptr;

  // write Lambda back to CTF
  std::vector<int64_t> lambdaIndices(UaiF->wrld->rank == 0 ? NvNo : 0);
  for (size_t i(0); i < lambdaIndices.size(); ++i) { lambdaIndices[i] = i; }
  lambdaF = new Tensor<>(1, &NvNo, Vbija->sym, *Vbija->wrld, "Lambda");
  lambdaF->write(lambdaIndices.size(), lambdaIndices.data(), lambdas.data());

  allocatedTensorArgument<>("SinglesHamiltonianEigenvalues", lambdaF);
}

void ThermalClusterDoublesAlgorithm::propagateAmplitudes(
  Tensor<real> &Sabij,
  const std::function<void(real, real &)> &propagator,
  Tensor<real> &S1FG
) {
  Transform<real, real> divBy(
    std::function<void(real, real &)>([](const real g, real &t){ t /= g; })
  );
  // half-open contraction weights for propagation
  divBy( (*ga)["a"], Sabij["abij"] );
  divBy( (*ga)["b"], Sabij["abij"] );
  divBy( (*gi)["i"], Sabij["abij"] );
  divBy( (*gi)["j"], Sabij["abij"] );
  // transform to singles-mode space
  Tensor<real> SFG(false, S1FG);
  SFG["FG"] = (*UaiF)["bjG"] * (*UaiF)["aiF"] * Sabij["abij"];
  // propagate
  Transform<real, real>(std::function<void(real, real &)>(propagator))(
    (*lambdaFG)["FG"], SFG["FG"]
  );
  // collect
  S1FG["FG"] -= SFG["FG"];
}

void ThermalClusterDoublesAlgorithm::diagonalizeDoublesAmplitudes() {
  LOG(1, getCapitalizedAbbreviation())
    << "diagonalizing doubles amplitudes at tau_1..." << std::endl;
  BlacsWorld world(Tabijn.back()->wrld->rank, Tabijn.back()->wrld->np);
  auto epsi(getTensorArgument<>("ThermalHoleEigenEnergies"));
  auto epsa(getTensorArgument<>("ThermalParticleEigenEnergies"));
  int NvNo(epsa->lens[0]*epsi->lens[0]);
  int scaTLens[2] = { NvNo, NvNo };
  auto scaTFG(NEW(ScaLapackMatrix<>, *Tabijn.back(), scaTLens, &world));

  // use ScaLapack routines to diagonalise the matrix U.Lambda.U^T
  auto scaU(NEW(ScaLapackMatrix<>, *scaTFG));
  ScaLapackHermitianEigenSystemDc<> eigenSystem(scaTFG, scaU);
  std::vector<real> lambdas(Tabijn.back()->lens[0]);
  eigenSystem.solve(lambdas.data());
  scaTFG = nullptr;

  if (isArgumentGiven("DoublesAmplitudesEigenvalues")) {
    // write Lambda back to CTF
    std::vector<int64_t> lambdaIndices(
      UaiF->wrld->rank == 0 ? lambdas.size() : 0
    );
    for (size_t i(0); i < lambdaIndices.size(); ++i) { lambdaIndices[i] = i; }
    auto lambdaL(
      new Tensor<>(
        1, scaTLens, Tabijn.front()->sym, *Tabijn.front()->wrld, "L"
      )
    );
    lambdaL->write(lambdaIndices.size(), lambdaIndices.data(), lambdas.data());
    allocatedTensorArgument<>("DoublesAmplitudesEigenvalues", lambdaF);
  }
}

