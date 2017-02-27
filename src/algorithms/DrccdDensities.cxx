#include <algorithms/DrccdDensities.hpp>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <util/CtfMachineTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace tcc;

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

  // call generic (dry and non-dry) run method
  run<CTF::Tensor<>, CtfMachineTensor<>>(epsi, epsa, ctfDabij, false);

  // convert hole occupancies from particle/hole operator N'_i = i^\dagger i
  // to electron N_i - c_i^\dagger c_i, where i^\dagger = c_i
  // TODO: implement below operation in tcc
  auto Ni(getTensorArgument<>("DrccdHoleOccupancies"));
  (*Ni)["i"] += 2.0;
}

void DrccdDensities::dryRun() {
  DryTensor<> *dryDabij(
    new DryTensor<>(
      *getTensorArgument<double, DryTensor<>>("DrccdDoublesAmplitudes")
    )
  );
  run<DryTensor<>, DryMachineTensor<>>(dryDabij, nullptr, nullptr, true);
}

template <typename T, typename MT>
void DrccdDensities::run(T *ctfEpsi, T *ctfEpsa, T *ctfDabij, const bool ) {
  auto machineTensorFactory(MT::Factory::create());
  auto tcc(Tcc<double>::create(machineTensorFactory));

  // Read the Drccd doubles amplitudes Tabij
  T *ctfTabij( getTensorArgument<double, T>("DrccdDoublesAmplitudes") );
  auto Tabij( tcc->createTensor(MT::create(*ctfTabij)) );

  // Read the required Coulomb integrals Vabij
  T *ctfVabij( getTensorArgument<double, T>("PPHHCoulombIntegrals") );
  auto Vabij( tcc->createTensor(MT::create(*ctfVabij)) );
  T *ctfVaibj( getTensorArgument<double, T>("PHPHCoulombIntegrals") );
  auto Vaibj( tcc->createTensor(MT::create(*ctfVaibj)) );
  T *ctfVijkl( getTensorArgument<double, T>("HHHHCoulombIntegrals") );
  auto Vijkl( tcc->createTensor(MT::create(*ctfVijkl)) );

  // convert energy denominators tensor into tcc tensor
  auto Dabij( tcc->createTensor(MT::create(*ctfDabij)) );
  delete ctfDabij;
  // same for the HF eigenvalues
  auto epsi( tcc->createTensor(MT::create(*ctfEpsi)) );
  auto epsa( tcc->createTensor(MT::create(*ctfEpsa)) );


  // create Lambda amplitudes and residuum of same shape as doubles amplitudes
  auto Labij( tcc->createTensor(Tabij, "Labij") );
  auto Rabij( tcc->createTensor(Tabij, "Rabij") );

  // compile Lambda equation iteration
  auto iterationOperation(
    tcc->compile( (
      (*Rabij)["abij"] <<= (*Vabij)["abij"],
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
      // note the sign from breaking up a hole line
      (*Dij)["ij"] <<= -2 * (*Tabij)["cdkj"] * (*Labij)["cdki"],
      (*Dab)["ab"] <<= +2 * (*Tabij)["cbkl"] * (*Labij)["cakl"],
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

  // build single particle energies (eigenvalues of HF - effective interaction)
  // this lets T + V_ne = Np*ep
  auto ei( tcc->createTensor(epsi, "ei") );
  auto ea( tcc->createTensor(epsa, "ea") );
  auto Veei( tcc->createTensor(epsi, "Veei") );
  auto Veea( tcc->createTensor(epsa, "Veea") );
  tcc->compile(
    (
      (*Veei)["i"] <<= 2 * (*Vijkl)["ijij"],
      (*Veei)["i"] -= (*Vijkl)["ijji"],
      (*Veea)["a"] <<= 2 * (*Vaibj)["ajaj"],
      (*Veea)["a"] -= (*Vabij)["aajj"],
      (*ei)["i"] <<= (*epsi)["i"],
      (*ei)["i"] -= (*Veei)["i"],
      (*ea)["a"] <<= (*epsa)["a"],
      (*ea)["a"] -= (*Veea)["a"]
    )
  )->execute();
  allocatedTensorArgument<double, T>(
    "DrccdOneBodyCoulombHoleEnergies",
    new T(Veei->template getMachineTensor<MT>()->tensor)
  );
  allocatedTensorArgument<double, T>(
    "DrccdOneBodyCoulombParticleEnergies",
    new T(Veea->template getMachineTensor<MT>()->tensor)
  );
  allocatedTensorArgument<double, T>(
    "DrccdOneBodyHoleEnergies",
    new T(ei->template getMachineTensor<MT>()->tensor)
  );
  allocatedTensorArgument<double, T>(
    "DrccdOneBodyParticleEnergies",
    new T(ea->template getMachineTensor<MT>()->tensor)
  );

  // build reduced two body density matrix
  // additionally, evaluate V_ee
  auto Gabij( tcc->createTensor(Vabij, "Gabij") );
  auto Vee( tcc->createTensor(std::vector<int>(), "Vee") );
  tcc-> compile(
    (
      // from Gamma^ab_ij
      (*Gabij)["abij"] <<= (*Tabij)["abij"],
      (*Gabij)["abij"] += 4*(*Tabij)["acik"]*(*Labij)["cdkl"]*(*Tabij)["dblj"],
      // from Gamma^aj_ib
      (*Gabij)["abij"] += 2 * (*Tabij)["acik"] * (*Labij)["bcjk"],
      // from Gamma^ib_aj
      (*Gabij)["abij"] += 2 * (*Tabij)["cbkj"] * (*Labij)["caki"],
      // from Gamma ^ij_ab
      (*Gabij)["abij"] += (*Labij)["abij"],
      // calcualte Coulomb energy beyond first order
      (*Vee)[""] <<= 2 * (*Gabij)["abij"] * (*Vabij)["abij"],
      // the last ineraction is also exchanged in drCCD: i.e. <0|1 V T|0>
      (*Vee)[""] -= (*Tabij)["abij"] * (*Vabij)["abji"]
    )
  )->execute();
  allocatedTensorArgument<double, T>(
    "DrccdTwoBodyPPHHDensiy",
    new T(Gabij->template getMachineTensor<MT>()->tensor)
  );
  allocatedTensorArgument<double, T>(
    "DrccdCoulombExpectationValue",
    new T(Vee->template getMachineTensor<MT>()->tensor)
  );
}

