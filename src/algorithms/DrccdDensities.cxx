#include <algorithms/DrccdDensities.hpp>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <util/CtfMachineTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace tcc;
using std::shared_ptr;

ALGORITHM_REGISTRAR_DEFINITION(DrccdDensities);

DrccdDensities::DrccdDensities(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

DrccdDensities::~DrccdDensities() {
}

void DrccdDensities::run() {
  // Read and bulld energy denominators
  CTF::Tensor<> *epsi( getTensorArgument<>("HoleEigenEnergies") );
  CTF::Tensor<> *epsa( getTensorArgument<>("ParticleEigenEnergies") );
  CTF::Tensor<> *ctfDabij(
    new CTF::Tensor<>(false, *getTensorArgument<>("DrccdDoublesAmplitudes"))
  );
  ctfDabij->set_name("Dabij");
  (*ctfDabij)["abij"] =  (*epsi)["i"];
  (*ctfDabij)["abij"] += (*epsi)["j"];
  (*ctfDabij)["abij"] -= (*epsa)["a"];
  (*ctfDabij)["abij"] -= (*epsa)["b"];
  CTF::Transform<double>(
    std::function<void(double &)>(
      [](double &d) { d = 1 / d; }
    )
  ) (
    (*ctfDabij)["abij"]
  );
  run<CTF::Tensor<>, CtfMachineTensor<>>(ctfDabij, false);
  auto Ni(getTensorArgument<>("DrccdHoleOccupancies"));
  (*Ni)["i"] -= 2.0;
  (*Ni)["i"] *= -1.0;
}

void DrccdDensities::dryRun() {
  DryTensor<> *dryDabij(
    new DryTensor<>(
      *getTensorArgument<double, DryTensor<>>("DrccdDoublesAmplitudes")
    )
  );
  run<DryTensor<>, DryMachineTensor<>>(dryDabij, true);
}

template <typename T, typename MT>
void DrccdDensities::run(T *ctfDabij, const bool dry) {
  auto machineTensorFactory(MT::Factory::create());
  auto tcc(Tcc<double>::create(machineTensorFactory));

  // Read the particle/hole Coulomb integrals Vabij
  T *ctfTabij( getTensorArgument<double, T>("DrccdDoublesAmplitudes") );
  auto Tabij( tcc->createTensor(MT::create(*ctfTabij)) );

  // Read the Drccd doubles amplitudes Tabij
  T *ctfVabij( getTensorArgument<double, T>("PPHHCoulombIntegrals") );
  auto Vabij( tcc->createTensor(MT::create(*ctfVabij)) );

  // convert energy denominators tensor into tcc tensor
  auto Dabij( tcc->createTensor(MT::create(*ctfDabij)) );
  delete ctfDabij;

  // create Lambda amplitudes and residuum of same shape as doubles amplitudes
  auto Labij( tcc->createTensor(Tabij, "Labij") );
  auto Rabij( tcc->createTensor(Tabij, "Rabij") );

  // compile Lambda equation iteration
  auto iterationOperation(
    tcc->compile( (
      (*Rabij)["abij"] <<= (*Vabij)["abij"],
      (*Rabij)["abij"] += (*Vabij)["abij"],
      (*Rabij)["abij"] += 2 * (*Labij)["acik"] * (*Vabij)["cbkj"],
      (*Rabij)["abij"] += 2 * (*Vabij)["acik"] * (*Labij)["cbkj"],
      (*Rabij)["abij"] += 4*(*Labij)["acik"]*(*Tabij)["cdkl"]*(*Vabij)["dblj"],
      (*Rabij)["abij"] += 4*(*Vabij)["acik"]*(*Tabij)["cdkl"]*(*Labij)["dblj"],
      (*Labij)["abij"] <<= (*Rabij)["abij"] * (*Dabij)["abij"]
    ) )
  );

  // execute iterations
  int maxIterations(getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS));
  for (int i(0); i < maxIterations; ++i) {
    LOG(0, "LambdaDRCCD") << "Iteration " << i << "..." << std::endl;
    iterationOperation->execute();
  }

  allocatedTensorArgument<double, T>(
    "DrccdLambdaDoublesAmplitudes",
    new T(Labij->template getMachineTensor<MT>()->tensor)
  );

  // build reduced one body density matrix
  // additionally, build number operator expectation values
  int Nv(Labij->lens[0]), No(Labij->lens[2]);
  auto Dij( tcc->createTensor(std::vector<int>({No,No}), "Dij") );
  auto Dab( tcc->createTensor(std::vector<int>({Nv,Nv}), "Dab") );
  auto Ni( tcc->createTensor(std::vector<int>({No}), "Ni") );
  auto Na( tcc->createTensor(std::vector<int>({Nv}), "Na") );
  tcc->compile(
    (
      (*Dij)["ij"] <<= 2 * (*Tabij)["cdkj"] * (*Labij)["cdki"],
      (*Dab)["ab"] <<= 2 * (*Tabij)["cbkl"] * (*Labij)["cakl"],
      (*Ni)["i"] <<= 2 * (*Dij)["ii"],
      (*Na)["a"] <<= 2 * (*Dab)["aa"]
    )
  )->execute();
  allocatedTensorArgument<double, T>(
    "DrccdOneBodyHHDensity", new T(Dij->template getMachineTensor<MT>()->tensor)
  );
  allocatedTensorArgument<double, T>(
    "DrccdOneBodyPPDensity", new T(Dab->template getMachineTensor<MT>()->tensor)
  );
  allocatedTensorArgument<double, T>(
    "DrccdHoleOccupancies", new T(Ni->template getMachineTensor<MT>()->tensor)
  );
  allocatedTensorArgument<double, T>(
    "DrccdParticleOccupancies",new T(Na->template getMachineTensor<MT>()->tensor)
  );
}

