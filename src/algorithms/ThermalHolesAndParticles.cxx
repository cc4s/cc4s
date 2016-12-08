// FIXME: bug at T with overlap but where not all occupied states are virtuals.
// i.e.: 0<=i<7 and 1<=a<30 in Mp2
#include <algorithms/ThermalHolesAndParticles.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

#include <algorithm>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(ThermalHolesAndParticles);

ThermalHolesAndParticles::ThermalHolesAndParticles(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ThermalHolesAndParticles::~ThermalHolesAndParticles() {
}

void ThermalHolesAndParticles::run() {
  orderStates();
  determineChemicalPotential();
  defineThermalHolesAndParticles();
  determineThermalOccupancies();

  delete epsp;
//  return 1.0 / (1.0 + std::exp(-(eps-mu)/kT));
  eigenStates.clear();
}


void ThermalHolesAndParticles::dryRun() {
}

void ThermalHolesAndParticles::orderStates() {
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  int N0(0), No(epsi->lens[0]), Nv(epsa->lens[0]), Np(No+Nv);

  // join (zero temperature) hole- and particle eigenenergies
  epsp = new Vector<>(Np, *epsi->wrld, "epsp");
  epsp->slice(&N0, &No, 0.0, *epsi, &N0, &No, 1.0);
  epsp->slice(&No, &Np, 0.0, *epsa, &N0, &Nv, 1.0);

  double *eigenEnergies(new double[Np]);
  epsp->read_all(eigenEnergies);

  // determine permutation of the states necessary to have ordered eigenvalues
  for (int p(0); p < Np; ++p) {
    eigenStates.push_back(std::pair<double, int>(eigenEnergies[p], p));
  }
  std::sort(eigenStates.begin(), eigenStates.end());
  int *sourceStates(new int[Np]);
  for (int p(0); p < Np; ++p) {
//    sourceStates[p] = eigenStates[p].second;
    sourceStates[eigenStates[p].second] = p;
  }
  // order the EigenEnergies vector
  // FIXME: why does CTF complain at the below call:
  // TF ERROR: please use other variant of permute function when the output is sparse
  // the output is not sparse...
//  epsp->permute(&sourceStates, 0.0, *epsp, 1.0);
  epsp->permute(0.0, *epsp, &sourceStates, 1.0);
  // check if they are really ordered
  epsp->read_all(eigenEnergies);
  for (int p(0); p < Np; ++p) {
    if (eigenEnergies[p] != eigenStates[p].first)
      throw new EXCEPTION("Failed to order EigenEnergies");
  }

  // order the states in the CoulombVertex accordingly
  Tensor<complex> *GammaGqr(getTensorArgument<complex>("CoulombVertex"));
  int *GammaPermutations[3] = { nullptr, sourceStates, sourceStates };
  LOG(1, "ThermalHolesAndParticles") << "Ordering states in CoulombVertex..." <<
    std::endl;
//  GammaGqr->permute(GammaPermutations, 0.0, *GammaGqr, 1.0);
  GammaGqr->permute(0.0, *GammaGqr, GammaPermutations, 1.0);
  // NOTE that this overwrites the CoulombVertex.
  // Maybe it's better to write ThermalCoulombVertex instead

  delete[] eigenEnergies;
  delete[] sourceStates;
}

void ThermalHolesAndParticles::determineChemicalPotential() {
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  int No(epsi->lens[0]);
  int Np(epsp->lens[0]);

  // get lower and upper bounds for the chemical potential
  double muLower(eigenStates[0].first);
  double muUpper(eigenStates[No].first);
  // actual number of electrons (states for spin restricted) in the system
  double NElectrons(0.5 * getRealArgument("Electrons"));
  kT = getRealArgument("Temperature");
  // current estimate for the chemical potential
  mu = 0.0;
  // expectation value of the number operator for current mu
  double N(0.0);
  double tolerance(getRealArgument("tolerance", 1e-8));
  int maxIterations(getIntegerArgument("maxIterations", 64));
  int iterations;
  for (iterations = 0; iterations < maxIterations; ++iterations) {
    mu = 0.5 * (muUpper + muLower);
    N = -NElectrons;
    for (int p(0); p < Np; ++p) {
      N += 1.0 / (1.0 + std::exp(+(eigenStates[p].first-mu)/kT));
    }
    if (N < -tolerance) {
      muLower = mu;
    } else if (N > +tolerance) {
      muUpper = mu;
    } else {
      break;
    }
  }
  if (iterations == maxIterations)
    throw new EXCEPTION("Failed to determine chemical potential.");
  LOG(1, "ThermalHolesAndParticles") << "Chemical potential=" << mu << std::endl;
  setRealArgument("ChemicalPotential", mu);
}

void ThermalHolesAndParticles::defineThermalHolesAndParticles() {
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  int No(epsi->lens[0]);
  int Np(epsp->lens[0]);

  // determine the thermal hole and particle index set
  double maxDelta(eigenStates[Np-1].first - eigenStates[0].first);
  double overlap(std::sqrt(32*maxDelta*kT));
  LOG(1, "ThermalHolesAndParticles") << "overlap of thermal particle and hole states=" << overlap << std::endl;
  int holeEnd(No), particleStart(-1);
  for (int p(0); p < Np; ++p) {
    if (eigenStates[p].first <= mu + overlap) holeEnd = p;
    if (eigenStates[p].first < mu - overlap) particleStart = p;
    LOG(4, "ThermalHolesAndParticles") << "eigenEnergies[" <<
      p << "<-" << eigenStates[p].second << "]=" <<
      eigenStates[p].first <<std::endl;
  }
  // holeEnd points to the last acceptable hole energy, we need to include it
  ++holeEnd;
  // particleStart points to the last unacceptable particle energy,
  // we need to exclude it
  ++particleStart;
  LOG(1, "ThermalHolesAndParticles") << "thermal hole states: 0<=i<" <<
    holeEnd << std::endl;
  LOG(1, "ThermalHolesAndParticles") << "thermal particle states: " <<
    particleStart << "<=a<" << Np << std::endl;

  Scalar<> chemicalPotential(mu, *epsp->wrld);

  // start and number of finite temperature virtual orbitals
  int ftN0(0), ftNv(Np-particleStart);
  // define thermal hole energies
  Vector<> *ftEpsi(new Vector<>(holeEnd, *epsp->wrld, "ftEpsi"));
  ftEpsi->slice(&ftN0, &holeEnd, 0.0, *epsp, &ftN0, &holeEnd, 1.0);
  (*ftEpsi)["i"] -= chemicalPotential[""];
  allocatedTensorArgument<>("ThermalHoleEigenEnergies", ftEpsi);

  // define thermal particle energies
  Vector<> *ftEpsa(new Vector<>(ftNv, *epsp->wrld, "ftEpsa"));
  ftEpsa->slice(&ftN0, &ftNv, 0.0, *epsp, &particleStart, &Np, 1.0);
  (*ftEpsa)["a"] -= chemicalPotential[""];
  allocatedTensorArgument<>("ThermalParticleEigenEnergies", ftEpsa);
}

void ThermalHolesAndParticles::determineThermalOccupancies() {
  // use positive temperatures for particles, negative ones for holes
  class occupancy {
  public:
    occupancy(double kT_): kT(kT_) { }
    void operator()(double &eps) {
      eps = 1.0 / (1.0 + std::exp(-eps/kT));
    }
  protected:
    double kT;
  };

  Tensor<> *ftEpsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> *ftNi(new Tensor<>(*ftEpsi));
  CTF::Transform<>( std::function<void(double &)>(occupancy(-kT)) ) (
    (*ftNi)["i"]
  );
  allocatedTensorArgument("ThermalHoleOccupancies", ftNi);

  Tensor<> *ftEpsa(getTensorArgument("ThermalParticleEigenEnergies"));
  Tensor<> *ftNa(new Tensor<>(*ftEpsa));
  CTF::Transform<>( std::function<void(double &)>(occupancy(+kT)) ) (
    (*ftNa)["a"]
  );
  allocatedTensorArgument("ThermalParticleOccupancies", ftNa);
}

