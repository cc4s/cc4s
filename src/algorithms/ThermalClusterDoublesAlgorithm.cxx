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

  // number of iterations for determining the amplitudes at each point in time
  int I( getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS) );
  // allocate doubles amplitudes on imaginary time grid
  TFn.resize(taus.size());
  TFGn.resize(taus.size());
  for (size_t n(0); n < taus.size(); ++n) {
    TFn[n] = NEW(CTF::Tensor<real>, false, *H0F);
    TFGn[n] = NEW(CTF::Tensor<real>, false, *VdFG);
  }
  std::vector<PTR(Tensor<real>)> SFn(taus.size()), SFGn(taus.size());
  for (size_t n(0); n < taus.size(); ++n) {
    SFn[n] = NEW(CTF::Tensor<real>, false, *H0F);
    SFGn[n] = NEW(CTF::Tensor<real>, false, *VdFG);
  }
  std::vector<real> energies(taus.size());

  // TODO: test then accept/dismiss renormalizations, for now disabled
  int R( getIntegerArgument("renormalizations", 0) );

  real energy;
  real accuracy(getRealArgument("accuracy", 1e-7));
  CTF::Tensor<real> T0F(false, *H0F);
  CTF::Tensor<real> T0FG(false, *VdFG);
  CTF::Tensor<real> S1F(false, *H0F);
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
        T0FG["FG"] = 0.0;
        S1FG["FG"] = 0.0;
        real direct(0.0), exchange(0.0), singles(0.0);
        for (size_t n(0); n <= N; ++n) {
          real tau1(scale*taus[n]);
          real DTau(tau1-tau0);
          // energy contribution from previously convolved amplitudes ST^I(tau_n-1)
          computeEnergyContribution(S1F, S1FG, DTau, direct, exchange, singles);

          // get T0FG at tau_n-1 and write previously convolved amplitudes
          // note T(0) is implicitly 0
          if (n > 0) {
            T0FG["FG"] = (*TFGn[n-1])["FG"];
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
          applyHamiltonian(T0F, T0FG, *TFn[n], *TFGn[n], DTau, S1F, S1FG);

          // energy contribution from convolved amplitudes ST^I(tau_n)
          computeEnergyContribution(S1F, S1FG, DTau, direct, exchange, singles);
          real a(S1FG.norm2());
          LOG(2, getCapitalizedAbbreviation())
            << "F_d=" << direct/tau1 << std::endl;
          LOG(2, getCapitalizedAbbreviation())
            << "F_x=" << exchange/tau1 << std::endl;
          if (getIntegerArgument("singlesEnergy", DEFAULT_SINGLES_ENERGY)) {
            LOG(2, getCapitalizedAbbreviation())
              << "F_s=" << singles/tau1 << std::endl;
          }
          LOG(2, getCapitalizedAbbreviation()) << "|T|=" << a << std::endl;
          tau0 = tau1;
        }
        (*SFGn[N])["FG"] = S1FG["FG"];
        energies[N] = energy = direct + exchange + singles;
        LOG(1, getCapitalizedAbbreviation()) << "F=" << energy/tau0 << std::endl;
        if (std::abs(1-lastEnergy/energy) < accuracy) break;
        lastEnergy = energy;
        real mixingRatio( getRealArgument("mixingRatio", 1.0) );
        for (size_t n(0); n <= N; ++n) {
          (*TFGn[n])["FG"] *= (1-mixingRatio);
          (*TFGn[n])["FG"] += mixingRatio * (*SFGn[n])["FG"];
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
  bool singlesEnergy(
    getIntegerArgument("singlesEnergy", DEFAULT_SINGLES_ENERGY)
  );

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
  // singles: one particle/hole pair F
  H0F = NEW(Tensor<real>, 1, std::vector<int>({NF}).data());
  // TODO: use Fock matrix for singles
//  (*H0F)["F"]  = (*UaiF)["ckF"] * (*ga)["c"] * (*gi)["k"]
//    * (*deltaai)["ck"] * (*epsa)["c"];

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

  real tdas(0.0);
  if (singlesEnergy) {
    // propagate singles
    Tensor<real> TF(*H0F);
    Transform<real, real>(std::function<void(real,real&)>(secondOrderIntegral))(
      (*lambdaF)["F"], TF["F"]
    );
    e[""] = -spins * TF["F"] * (*H0F)["F"];
    tdas = e.get_val();
    LOG(1, getCapitalizedAbbreviation()) << "TDA F_s=" << tdas << std::endl;
  }

  if (isArgumentGiven("SinglesHamiltonianWeights")) {
    auto singlesWeights(new Tensor<real>(1, std::vector<int>({NF}).data()) );
    (*singlesWeights)["F"] = 0.5 * spins*spins * TFG["FF"] * (*VdFG)["FF"];
    allocatedTensorArgument<real>("SinglesHamiltonianWeights", singlesWeights);
  }

  return tdad+tdax+tdas;
}

void ThermalClusterDoublesAlgorithm::computeEnergyContribution(
  Tensor<real> &T1F, Tensor<real> &T2FG, const real DTau,
  real &direct, real &exchange, real &singles
) {
  real spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  Scalar<real> energy;
  // combine singles and doubles
  Tensor<real> TFG(T2FG);
  bool singlesEnergy(
    getIntegerArgument("singlesEnergy", DEFAULT_SINGLES_ENERGY)
  );
  if (singlesEnergy) {
    TFG["FG"] += T1F["F"] * T1F["G"];
  }

  // direct term
  energy[""] = 0.5 * spins * spins * DTau/2 * TFG["FG"] * (*VdFG)["FG"];
  direct += energy.get_val();
  // exchange term
  energy[""] = (-0.5) * spins * DTau/2 * TFG["FG"] * (*VxFG)["FG"];
  exchange += energy.get_val();
  if (singlesEnergy) {
    // singles term
    energy[""] = spins * DTau/2 * T1F["F"] * (*H0F)["F"];
    singles += energy.get_val();
  }
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
  // unperturbed propatation is diagonal
  (*Hbjai)["bjbj"] += (*epsa)["b"];
  (*Hbjai)["bjbj"] -= (*epsi)["j"];
  // bubble from H_1 has contraction weights
  (*Hbjai)["bjai"] = spins * (*Vbija)["bija"];
  (*Hbjai)["bjai"] *= (*ga)["a"];
  (*Hbjai)["bjai"] *= (*ga)["b"];
  (*Hbjai)["bjai"] *= (*gi)["i"];
  (*Hbjai)["bjai"] *= (*gi)["j"];

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

void ThermalClusterDoublesAlgorithm::diagonalizeDoublesAmplitudes() {
  LOG(1, getCapitalizedAbbreviation())
    << "diagonalizing doubles amplitudes at tau_1..." << std::endl;
  BlacsWorld world(TFGn.back()->wrld->rank, TFGn.back()->wrld->np);
  auto scaTFG(NEW(ScaLapackMatrix<>, *TFGn.back(), TFGn.back()->lens, &world));

  // use ScaLapack routines to diagonalise the matrix U.Lambda.U^T
  auto scaU(NEW(ScaLapackMatrix<>, *scaTFG));
  ScaLapackHermitianEigenSystemDc<> eigenSystem(scaTFG, scaU);
  std::vector<real> lambdas(TFGn.back()->lens[0]);
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
        1, TFGn.front()->lens, TFGn.front()->sym, *TFGn.front()->wrld, "L"
      )
    );
    lambdaL->write(lambdaIndices.size(), lambdaIndices.data(), lambdas.data());
    allocatedTensorArgument<>("DoublesAmplitudesEigenvalues", lambdaF);
  }
}

