#include <algorithms/ClusterDoublesAlgorithm.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <tcc/Tcc.hpp>
#include <util/CtfMachineTensor.hpp>
#include <Options.hpp>

using namespace CTF;
using namespace cc4s;

ClusterDoublesAlgorithm::ClusterDoublesAlgorithm(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ClusterDoublesAlgorithm::~ClusterDoublesAlgorithm() {
}

void ClusterDoublesAlgorithm::run() {
  Data *Vabij(getArgumentData("PPHHCoulombIntegrals"));
  TensorData<double> *realVabij(dynamic_cast<TensorData<double> *>(Vabij));
  double e(0.0);
  if (realVabij) {
    e = run<double>();
  } else {
    e = std::real( run<complex>() );
  }

  std::stringstream energyName;
  energyName << getAbbreviation() << "Energy";
  setRealArgument(energyName.str(), e);
}

void ClusterDoublesAlgorithm::dryRun() {
  Data *Vabij(getArgumentData("PPHHCoulombIntegrals"));
  TensorData<double> *realVabij(dynamic_cast<TensorData<double> *>(Vabij));
  double e(0.0);
  if (realVabij) {
    e = dryRun<double>();
  } else {
    e = std::real( dryRun<complex>() );
  }

  std::stringstream energyName;
  energyName << getAbbreviation() << "Energy";
  setRealArgument(energyName.str(), e);
}


template <typename F>
F ClusterDoublesAlgorithm::run() {
  Mixer<F> *TabijMixer(createDoublesMixer<F>());

  // Iteration for determining the doubles amplitudes Tabij
  // and the energy e
  int maxIterationsCount(
    getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS)
  );
  F e(0);
  for (int i(0); i < maxIterationsCount; ++i) {
    LOG(0, getCapitalizedAbbreviation()) << "iteration: " << i+1 << std::endl;
    // call the iterate of the actual algorithm, which is still left open here
    // no singles amplitudes used in ClusterDoublesAlgorithm, such as CCD.
    iterate(i, nullptr, TabijMixer);
    e = calculateEnergy(TabijMixer);
  }

  // FIXME: support energy calculation from given amplitdues.
  return e;
}


template <typename F>
F ClusterDoublesAlgorithm::dryRun() {
  // Read the Coulomb Integrals Vabij required for the energy
  getTensorArgument<F, DryTensor<F>>("PPHHCoulombIntegrals");

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<F> *epsi(
    getTensorArgument<F, DryTensor<F>>("HoleEigenEnergies")
  );
  DryTensor<F> *epsa(
    getTensorArgument<F, DryTensor<F>>("ParticleEigenEnergies")
  );

  std::stringstream amplitudesName;
  amplitudesName << getAbbreviation() << "DoublesAmplitudes";

  // Instantiate mixer for the doubles amplitudes, by default use the linear one
  std::string mixerName(getTextArgument("mixer", "LinearMixer"));
  Mixer<F> *TabijMixer = MixerFactory<F>::create(mixerName, this);
  if (!TabijMixer) {
    std::stringstream stringStream;
    stringStream << "Mixer not implemented: " << mixerName;
    throw new EXCEPTION(stringStream.str());
  }
  // TODO: implement DryTensor in mixers
  if (mixerName != "LinearMixer") {
    LOG(0, getCapitalizedAbbreviation())
      << "Warning: dry run not implemented for " << mixerName
      << ", assuming the same memory usage." << std::endl;
  }

  DryTensor<F> *Tabij;
  {
    // Allocate the doubles amplitudes and append it to the mixer
    int No(epsi->lens[0]);
    int Nv(epsa->lens[0]);
    int syms[] = { NS, NS, NS, NS };
    int vvoo[] = { Nv, Nv, No, No };
    Tabij = new DryTensor<F>(4, vvoo, syms, SOURCE_LOCATION);
    allocatedTensorArgument<F>(amplitudesName.str(), Tabij);
    // FIXME: no dry mixer exists yet
//    TabijMixer->append(Tabij);
  }

  // Allocate the energy e
  DryScalar<F> energy(SOURCE_LOCATION);

  getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS);

  // Call the dry iterate of the actual algorithm, which is left open here
  // pass the doubles amplitudes but no singles amplitudes
  dryIterate(nullptr, Tabij);

  return 0;
}


// default dryIterate
void ClusterDoublesAlgorithm::dryIterate(
  DryTensor<double> *Tai, DryTensor<double> *Tabij
) {
  LOG(0, "CluserDoubles") << "Dry run for iteration not given for "
    << getAbbreviation() << std::endl;
}

void ClusterDoublesAlgorithm::dryIterate(
  DryTensor<complex> *Tai, DryTensor<complex> *Tabij
) {
  LOG(0, "CluserDoubles") << "Dry run for iteration not given for "
    << getAbbreviation() << std::endl;
}


template <typename F>
Mixer<F> *ClusterDoublesAlgorithm::createDoublesMixer() {
  // Instantiate mixer for the doubles amplitudes, by default use the linear one
  std::string mixerName(getTextArgument("mixer", "LinearMixer"));
  Mixer<F> *TabijMixer( MixerFactory<F>::create(mixerName, this));
  if (!TabijMixer) {
    std::stringstream stringStream;
    stringStream << "Mixer not implemented: " << mixerName;
    throw new EXCEPTION(stringStream.str());
  }

  if (isArgumentGiven("startingDoublesAmplitudes")) {
    Tensor<F> *Tabij(getTensorArgument<F>("startingDoublesAmplitudes"));
    // use given amplitudes as initial amplitudes
    TabijMixer->append(*Tabij);
  } else {
    // Read the Coulomb Integrals Vabij
    Tensor<F> *Vabij(getTensorArgument<F>("PPHHCoulombIntegrals"));
    // copy its shape only to get zeros as initial amplitudes
    Tensor<F> Tabij(false, *Vabij);
    Tabij.set_name("Tabij");
    TabijMixer->append(Tabij);
  }
  // The amplitudes will from now on be managed by the mixer
  return TabijMixer;
}


template <typename F>
void ClusterDoublesAlgorithm::storeDoublesAmplitudes(Mixer<F> *TabijMixer) {
  if (isArgumentGiven(getDoublesAmplitudesName())) {
    allocatedTensorArgument<F>(
      getDoublesAmplitudesName(), new Tensor<F>(TabijMixer->getNext())
    );
  }
}

template <typename F>
F ClusterDoublesAlgorithm::calculateEnergy(Mixer<F> *TabijMixer) {
  // get the Coulomb integrals to compute the energy
  Tensor<F> *Vijab(getTensorArgument<F>("HHPPCoulombIntegrals"));

  // allocate energy
  Scalar<F> energy(*Vijab->wrld);
  energy.set_name("energy");

  // get amplitudes from the mixer
  Tensor<F> *Tabij(&TabijMixer->getNext());
  Tabij->set_name("Tabij");

  // direct term
  energy[""] = +2.0 * (*Tabij)["abij"] * (*Vijab)["ijab"];
  F dire(energy.get_val());
  // exchange term
  energy[""] = -1.0 * (*Tabij)["abij"] * (*Vijab)["ijba"];
  F exce(energy.get_val());
  F e(dire + exce);
  LOG(0, getCapitalizedAbbreviation()) << "e=" << e << std::endl;
  LOG(1, getCapitalizedAbbreviation()) << "dir=" << dire << std::endl;
  LOG(1, getCapitalizedAbbreviation()) << "exc=" << exce << std::endl;
  return e;
}

template <typename F>
void ClusterDoublesAlgorithm::doublesAmplitudesFromResiduum(
  CTF::Tensor<F> &Rabij
) {
  Tensor<F> *Vabij(getTensorArgument<F>("PPHHCoulombIntegrals"));
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));

  // convert to type F (either complex or double)
  Tensor<F> Fepsi(1, &epsi->lens[0], epsi->sym, *epsi->wrld, "Fepsi");
  // NOTE: just copies if both arguments are real
  toComplexTensor(*epsi, Fepsi);
  Tensor<F> Fepsa(1, &epsa->lens[0], epsa->sym, *epsa->wrld, "Fepsa");
  toComplexTensor(*epsa, Fepsa);

  // create excitation energy tensor
  Tensor<F> Dabij(false, *Vabij);
  Dabij.set_name("Dabij");
  Dabij["abij"] =  Fepsi["i"];
  Dabij["abij"] += Fepsi["j"];
  Dabij["abij"] -= Fepsa["a"];
  Dabij["abij"] -= Fepsa["b"];

  // TODO:
  // levelshifting can be implemented here

  // use transform to divide Tabij by Dabij
  CTF::Transform<F, F>(
    std::function<void(F, F &)>(
      [](F dabij, F &t) {
        t = t / dabij;
      }
    )
  ) (
    Dabij["abij"], Rabij["abij"]
  );
}

void ClusterDoublesAlgorithm::dryDoublesAmplitudesFromResiduum(
  cc4s::DryTensor<> &Rabij
) {
  // Build Dabij
  DryTensor<> Dabij(Rabij, SOURCE_LOCATION);
}

// instantiate
template
void ClusterDoublesAlgorithm::doublesAmplitudesFromResiduum(
  CTF::Tensor<double> &Rabij
);

template
void ClusterDoublesAlgorithm::doublesAmplitudesFromResiduum(
  CTF::Tensor<complex> &Rabij
);


std::string ClusterDoublesAlgorithm::getCapitalizedAbbreviation() {
  std::string capitalizedAbbreviation(getAbbreviation());
  std::transform(
    capitalizedAbbreviation.begin(), capitalizedAbbreviation.end(),
    capitalizedAbbreviation.begin(), ::toupper
  );
  return capitalizedAbbreviation;
}

std::string ClusterDoublesAlgorithm::getDoublesAmplitudesName() {
  std::stringstream amplitudesName;
  amplitudesName << getAbbreviation() << "DoublesAmplitudes";
  return amplitudesName.str();
}

