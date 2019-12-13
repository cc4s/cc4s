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
  if (getArgumentData("ChemicalPotential")->getStage() != Data::Stage::MENTIONED) {
    determineNumberOfElectrons();
  } else {
    determineChemicalPotential();
  }
  defineThermalHolesAndParticles();
  determineThermalOccupancies();
  if (isArgumentGiven("ThermalParticleHoleOverlap")) {
    computeParticleHoleOverlap();
  }

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
  int Np(epsp->lens[0]);

  // get lower and upper bounds for the chemical potential
  double muLower(2*eigenStates[0].first-eigenStates[Np-1].first);
  double muUpper(2*eigenStates[Np-1].first-eigenStates[0].first);
  // actual number of electrons (states for spin restricted) in the system
  double spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  double NElectrons(getRealArgument("Electrons") / spins);
  kT = getRealArgument("Temperature");
  LOG(1,"FT")
    << "searching for mu in [" << muLower << "," << muUpper << "]" << std::endl;
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
  LOG(1, "ThermalHolesAndParticles") << "Chemical potential=" << mu << std::endl;
  if (iterations == maxIterations)
    throw new EXCEPTION("Failed to determine chemical potential.");
  setRealArgument("ChemicalPotential", mu);
}

void ThermalHolesAndParticles::determineNumberOfElectrons() {
  int Np(epsp->lens[0]);
  mu = getRealArgument("ChemicalPotential");
  // actual number of electrons (states for spin restricted) in the system
  double spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);
  kT = getRealArgument("Temperature");
  // expectation value of the number operator for current mu
  double N(0.0);
  for (int p(0); p < Np; ++p) {
    N += 1.0 / (1.0 + std::exp(+(eigenStates[p].first-mu)/kT));
  }
  N *= spins;
  LOG(1, "ThermalHolesAndParticles") << "Number of electrons according to reference=" << N << std::endl;
  setRealArgument("Electrons", N);
}

void ThermalHolesAndParticles::defineThermalHolesAndParticles() {
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  int No(epsi->lens[0]);
  int Np(epsp->lens[0]);

  // determine the thermal hole and particle index set
  double maxDelta(
    getRealArgument(
      "maxDenominator", eigenStates[Np-1].first - eigenStates[0].first
    )
  );
  double overlap(std::sqrt(32*maxDelta*kT));
  double minOccupancy( getRealArgument("minOccupancy", 1e-5) );
  int holeEnd(No), particleStart(-1);
  LOG(1, "ThermalHolesAndParticles") <<
    "max overlap of thermal particle and hole states=" << overlap << std::endl;
  for (int p(0); p < Np; ++p) {
    if (
      (eigenStates[p].first <= mu + overlap) &&
      (1.0 / (1.0 + std::exp(+(eigenStates[p].first-mu)/kT)) >= minOccupancy)
    ) holeEnd = p;
    if (
      (eigenStates[p].first < mu - overlap) ||
      (1.0 / (1.0 + std::exp(-(eigenStates[p].first-mu)/kT)) < minOccupancy)
    ) particleStart = p;
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
  allocatedTensorArgument<>("ThermalHoleEigenEnergies", ftEpsi);

  // define thermal particle energies
  Vector<> *ftEpsa(new Vector<>(ftNv, *epsp->wrld, "ftEpsa"));
  ftEpsa->slice(&ftN0, &ftNv, 0.0, *epsp, &particleStart, &Np, 1.0);
  allocatedTensorArgument<>("ThermalParticleEigenEnergies", ftEpsa);
}

void ThermalHolesAndParticles::determineThermalOccupancies() {
  // use positive temperatures for particles, negative ones for holes
  double minOccupancy( getRealArgument("minOccupancy", 1e-5) );
  class occupancy {
  public:
    occupancy(double kT_, double minf_): kT(kT_), minf(minf_) { }
    void operator()(double &eps) {
      double f(1.0 / (1.0 + std::exp(-eps/kT)));
      double invf(1.0 / (1.0 + std::exp(+eps/kT)));
      if (f < minf) {
        eps = 0.0;
      } else if (invf < minf) {
        eps = 1.0;
      } else {
        eps = f;
      }
    }
  protected:
    double kT, minf;
  };

  Tensor<> *ftEpsi(getTensorArgument("ThermalHoleEigenEnergies"));
  Tensor<> *ftNi(new Tensor<>(*ftEpsi));
  CTF::Transform<>(std::function<void(double &)>(occupancy(-kT,minOccupancy))) (
    (*ftNi)["i"]
  );
  allocatedTensorArgument("ThermalHoleOccupancies", ftNi);

  Tensor<> *ftEpsa(getTensorArgument("ThermalParticleEigenEnergies"));
  Tensor<> *ftNa(new Tensor<>(*ftEpsa));
  CTF::Transform<>(std::function<void(double &)>(occupancy(+kT,minOccupancy))) (
    (*ftNa)["a"]
  );
  allocatedTensorArgument("ThermalParticleOccupancies", ftNa);
}

void ThermalHolesAndParticles::computeParticleHoleOverlap() {
  auto epsi( getTensorArgument<>("ThermalHoleEigenEnergies") );
  int No(epsi->lens[0]);
  auto epsa( getTensorArgument<>("ThermalParticleEigenEnergies") );
  int Nv(epsa->lens[0]);
  Tensor<real> *deltaai(new Tensor<real>(2, std::vector<int>({Nv,No}).data()));
  auto GammaFqr( getTensorArgument<complex>("CoulombVertex") );
  int Np(GammaFqr->lens[1]);
  // there as many entries in delta^a_i as there is overlap between No and Nv
  size_t elementsCount(epsi->wrld->rank == 0 ? No+Nv-Np : 0);
  std::vector<real> elements(elementsCount);
  std::vector<int64_t> indices(elementsCount);
  for (size_t a(0); a < elementsCount; ++a) {
    elements[a] = 1.0;
    // delta^a_i = 1.0 iff a and i refer to the same index q
    int i = No - elementsCount + a;
    // assign the appropriate position
    indices[a] = a + Nv*i;
  }
  deltaai->write(elementsCount, indices.data(), elements.data());
  allocatedTensorArgument("ThermalParticleHoleOverlap", deltaai);
}

