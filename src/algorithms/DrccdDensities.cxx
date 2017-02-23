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

  // create Lambda operator amplitudes of same shape as doubles amplitudes
  auto Labij( tcc->createTensor(Tabij, "Labij") );

  // compile Lambda equation iteration
  auto iterationOperation(
    tcc->compile( (
      (*Labij)["abij"] <<= (*Vabij)["abij"],
      (*Labij)["abij"] += (*Vabij)["abij"],
      (*Labij)["abij"] += 2 * (*Labij)["acik"] * (*Vabij)["cbkj"],
      (*Labij)["abij"] += 2 * (*Vabij)["acik"] * (*Labij)["cbkj"],
      (*Labij)["abij"] += 4*(*Labij)["acik"]*(*Tabij)["cdkl"]*(*Vabij)["dblj"],
      (*Labij)["abij"] += 4*(*Vabij)["acik"]*(*Tabij)["cdkl"]*(*Labij)["dblj"],
      (*Labij)["abij"] <<= (*Labij)["abij"] * (*Dabij)["abij"]
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
}

