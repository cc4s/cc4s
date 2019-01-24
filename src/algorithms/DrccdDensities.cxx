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
  run<CTF::Tensor<real>, CtfMachineTensor<real>>(epsi, epsa, ctfDabij, false);

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
  run<DryTensor<real>, DryMachineTensor<real>>(dryDabij, nullptr, nullptr, true);
}

template <typename T, typename MT>
void DrccdDensities::run(T *ctfEpsi, T *ctfEpsa, T *ctfDabij, const bool ) {
  typedef typename MT::TensorEngine Engine;
  typedef tcc::Tcc<Engine> TCC;

  // Read the Drccd doubles amplitudes Tabij
  T *ctfTabij( getTensorArgument<double, T>("DrccdDoublesAmplitudes") );
  auto Tabij( tcc::Tensor<real,Engine>::create(*ctfTabij) );

  // Read the required Coulomb integrals Vabij
  T *ctfVabij( getTensorArgument<double, T>("PPHHCoulombIntegrals") );
  auto Vabij( tcc::Tensor<real,Engine>::create(*ctfVabij) );
  T *ctfVaibj( getTensorArgument<double, T>("PHPHCoulombIntegrals") );
  auto Vaibj( tcc::Tensor<real,Engine>::create(*ctfVaibj) );
  T *ctfVijkl( getTensorArgument<double, T>("HHHHCoulombIntegrals") );
  auto Vijkl( tcc::Tensor<real,Engine>::create(*ctfVijkl) );

  // convert energy denominators tensor into tcc tensor
  auto Dabij( tcc::Tensor<real,Engine>::create(*ctfDabij) );
  delete ctfDabij;
  // same for the HF eigenvalues
  auto epsi( tcc::Tensor<real,Engine>::create(*ctfEpsi) );
  auto epsa( tcc::Tensor<real,Engine>::create(*ctfEpsa) );


  // create Lambda amplitudes and residuum of same shape as doubles amplitudes
  auto Labij( TCC::template tensor(Tabij, "Labij") );
  auto Rabij( TCC::template tensor(Tabij, "Rabij") );

  double spins(getIntegerArgument("unrestricted", 0) ? 1.0 : 2.0);

  // compile Lambda equation iteration
  auto iterationOperation(
    (
      (*Rabij)["abij"] <<= (*Vabij)["abij"],
      (*Rabij)["abij"] += spins * (*Labij)["acik"] * (*Vabij)["cbkj"],
      (*Rabij)["abij"] += spins * (*Vabij)["acik"] * (*Labij)["cbkj"],
      (*Rabij)["abij"] += spins*spins*(*Labij)["acik"]*(*Tabij)["cdkl"]*(*Vabij)["dblj"],
      (*Rabij)["abij"] += spins*spins*(*Vabij)["acik"]*(*Tabij)["cdkl"]*(*Labij)["dblj"],
      (*Labij)["abij"] <<= (*Rabij)["abij"] * (*Dabij)["abij"]
    )->compile()
  );

  // execute iterations
  int maxIterations(getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS));
  for (int i(0); i < maxIterations; ++i) {
    LOG(0, "LambdaDRCCD") << "Iteration " << i << "..." << std::endl;
    iterationOperation->execute();
  }

  allocatedTensorArgument<double, T>(
    "DrccdLambdaDoublesAmplitudes",
    new T(Labij->getMachineTensor()->tensor)
  );

  // build reduced one body density matrix
  // additionally, build number operator expectation values
  auto Nv(Labij->lens[0]), No(Labij->lens[2]);
  auto Dij( TCC::template tensor(std::vector<size_t>({No,No}), "Dij") );
  auto Dab( TCC::template tensor(std::vector<size_t>({Nv,Nv}), "Dab") );
  auto Ni( TCC::template tensor(std::vector<size_t>({No}), "Ni") );
  auto Na( TCC::template tensor(std::vector<size_t>({Nv}), "Na") );
  (
    // note the sign from breaking up a hole line
    (*Dij)["ij"] <<= -spins * (*Tabij)["cdkj"] * (*Labij)["cdki"],
    (*Dab)["ab"] <<= +spins * (*Tabij)["cbkl"] * (*Labij)["cakl"],
    (*Ni)["i"] <<= spins * (*Dij)["ii"],
    (*Na)["a"] <<= spins * (*Dab)["aa"]
  )->compile()->execute();
  allocatedTensorArgument<double, T>(
    "DrccdOneBodyHHDensity", new T(Dij->getMachineTensor()->tensor)
  );
  allocatedTensorArgument<double, T>(
    "DrccdOneBodyPPDensity", new T(Dab->getMachineTensor()->tensor)
  );
  allocatedTensorArgument<double, T>(
    "DrccdHoleOccupancies", new T(Ni->getMachineTensor()->tensor)
  );
  allocatedTensorArgument<double, T>(
    "DrccdParticleOccupancies",new T(Na->getMachineTensor()->tensor)
  );

  // build single particle energies (eigenvalues of HF - effective interaction)
  // this lets T + V_ne = Np*ep
  auto ei( TCC::template tensor(epsi, "ei") );
  auto ea( TCC::template tensor(epsa, "ea") );
  auto Veei( TCC::template tensor(epsi, "Veei") );
  auto Veea( TCC::template tensor(epsa, "Veea") );
  (
    (*Veei)["i"] <<= 0.5*spins*spins * (*Vijkl)["ijij"],
    (*Veei)["i"] -= 0.5*spins * (*Vijkl)["ijji"],
    (*Veea)["a"] <<= 0.5*spins*spins * (*Vaibj)["ajaj"],
    (*Veea)["a"] -= 0.5*spins * (*Vabij)["aajj"],
    (*ei)["i"] <<= (*epsi)["i"],
    (*ei)["i"] -= (*Veei)["i"],
    (*ea)["a"] <<= (*epsa)["a"],
    (*ea)["a"] -= (*Veea)["a"]
  )->compile()->execute();
  allocatedTensorArgument<double, T>(
    "DrccdOneBodyCoulombHoleEnergies",
    new T(Veei->getMachineTensor()->tensor)
  );
  allocatedTensorArgument<double, T>(
    "DrccdOneBodyCoulombParticleEnergies",
    new T(Veea->getMachineTensor()->tensor)
  );
  allocatedTensorArgument<double, T>(
    "DrccdOneBodyHoleEnergies",
    new T(ei->getMachineTensor()->tensor)
  );
  allocatedTensorArgument<double, T>(
    "DrccdOneBodyParticleEnergies",
    new T(ea->getMachineTensor()->tensor)
  );

  // build reduced two body density matrix
  // additionally, evaluate V_ee
  auto Gabij( TCC::template tensor(Vabij, "Gabij") );
  auto Vee( TCC::template tensor(std::vector<size_t>(), "Vee") );
  (
    // from Gamma^ab_ij
    (*Gabij)["abij"] <<= (*Tabij)["abij"],
    (*Gabij)["abij"] += spins*spins * (*Tabij)["acik"]*(*Labij)["cdkl"]*(*Tabij)["dblj"],
    // from Gamma^aj_ib
    (*Gabij)["abij"] += spins * (*Tabij)["acik"] * (*Labij)["bcjk"],
    // from Gamma^ib_aj
    (*Gabij)["abij"] += spins * (*Tabij)["cbkj"] * (*Labij)["caki"],
    // from Gamma ^ij_ab
    (*Gabij)["abij"] += (*Labij)["abij"],
    // calcualte Coulomb energy beyond first order
    (*Vee)[""] <<= 0.5*spins*spins * (*Gabij)["abij"] * (*Vabij)["abij"],
    // the last ineraction is also exchanged in drCCD: i.e. <0|1 V T|0>
    (*Vee)[""] -= 0.5*spins * (*Tabij)["abij"] * (*Vabij)["abji"]
  )->compile()->execute();
  allocatedTensorArgument<double, T>(
    "DrccdTwoBodyPPHHDensiy",
    new T(Gabij->getMachineTensor()->tensor)
  );
  allocatedTensorArgument<double, T>(
    "DrccdCoulombExpectationValue",
    new T(Vee->getMachineTensor()->tensor)
  );
}

