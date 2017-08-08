#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <Options.hpp>

using namespace CTF;
using namespace cc4s;

ClusterSinglesDoublesAlgorithm::ClusterSinglesDoublesAlgorithm(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ClusterSinglesDoublesAlgorithm::~ClusterSinglesDoublesAlgorithm() {
}


void ClusterSinglesDoublesAlgorithm::run() {
  Data *Vabij(getArgumentData("PPHHCoulombIntegrals"));
  TensorData<double> *realVabij(dynamic_cast<TensorData<double> *>(Vabij));
  double e(0.0);
  if (realVabij) {
    e = run<double>();
  } else {
    e = std::real( run<complex>() );
  }
  setRealArgument(getDataName("", "Energy"), e);
}


template <typename F>
F ClusterSinglesDoublesAlgorithm::run() {
  int Nv(getTensorArgument<>("ParticleEigenEnergies")->lens[0]);
  int No(getTensorArgument<>("HoleEigenEnergies")->lens[0]);
  Mixer<F> *TaiMixer(
    createMixer<F>("Singles", std::vector<int>{{Nv,No}})
  );
  Mixer<F> *TabijMixer(
    createMixer<F>("Doubles", std::vector<int>{{Nv,Nv,No,No}})
  );

  // Iteration for determining the amplitudes
  int maxIterationsCount(
    getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS)
  );

  F e(0);
  for (int i(0); i < maxIterationsCount; ++i) {
    LOG(0, getCapitalizedAbbreviation()) << "iteration: " << i+1 << std::endl;
    // call the iterate of the actual algorithm, which is still left open here
    iterate(i, TaiMixer, TabijMixer);

    e = calculateEnergy(TaiMixer, TabijMixer);
  }

  if (maxIterationsCount == 0) {
    LOG(0, getCapitalizedAbbreviation()) <<
      "computing energy from given amplitudes" << std::endl;
    e = calculateEnergy(TaiMixer, TabijMixer);
  }

  storeAmplitudes(TaiMixer, "Singles");
  storeAmplitudes(TabijMixer, "Doubles");

  return e;
}



template <typename F>
F ClusterSinglesDoublesAlgorithm::calculateEnergy(
  Mixer<F> *TaiMixer, Mixer<F> *TabijMixer
) {
  // get the Coulomb integrals to compute the energy
  Tensor<F> *Vijab(getTensorArgument<F>("HHPPCoulombIntegrals"));

  // allocate energy
  Scalar<F> energy(*Vijab->wrld);
  energy.set_name("energy");

  // singles amplitudes are optional
  Tensor<F> *Tai( TaiMixer ? &TaiMixer->getNext() : nullptr );
  Tensor<F> *Tabij(&TabijMixer->getNext());

  // direct term
  energy[""] =  +2.0 * (*Tabij)["abij"] * (*Vijab)["ijab"];
  if (Tai) {
    energy[""] += +2.0 * (*Tai)["ai"] * (*Tai)["bj"] * (*Vijab)["ijab"];
  }
  F dire(energy.get_val());
  // exchange term
  energy[""] =  -1.0 * (*Tabij)["abij"] * (*Vijab)["ijba"];
  if (Tai) {
    energy[""] += -1.0 * (*Tai)["ai"] * (*Tai)["bj"] * (*Vijab)["ijba"];
  }
  F exce(energy.get_val());
  F e(dire + exce);
  LOG(0, getCapitalizedAbbreviation()) << "e=" << e << std::endl;
  LOG(1, getCapitalizedAbbreviation()) << "dir=" << dire << std::endl;
  LOG(1, getCapitalizedAbbreviation()) << "exc=" << exce << std::endl;
  return e;
}


template <typename F>
Mixer<F> *ClusterSinglesDoublesAlgorithm::createMixer(
  const std::string &type, std::vector<int> shape
) {
  // Instantiate mixer for the doubles amplitudes, by default use the linear one
  std::string mixerName(getTextArgument("mixer", "LinearMixer"));
  Mixer<F> *mixer( MixerFactory<F>::create(mixerName, this));
  if (!mixer) {
    std::stringstream stringStream;
    stringStream << "Mixer not implemented: " << mixerName;
    throw new EXCEPTION(stringStream.str());
  }

  std::stringstream initialDataName;
  initialDataName << "initial" << type << "Amplitudes";
  if (isArgumentGiven(initialDataName.str())) {
    Tensor<F> *T(getTensorArgument<F>(initialDataName.str()));
    // use given amplitudes as initial amplitudes
    mixer->append(*T);
  } else {
    Tensor<F> *V(getTensorArgument<F>("PPHHCoulombIntegrals"));
    std::vector<int> syms(shape.size(), NS);
    Tensor<F> T(shape.size(), shape.data(), syms.data(), *V->wrld, "T");
    mixer->append(T);
  }
  // The amplitudes will from now on be managed by the mixer
  return mixer;
}


template <typename F>
void ClusterSinglesDoublesAlgorithm::storeAmplitudes(
  Mixer<F> *mixer, const std::string &type
) {
  if (isArgumentGiven(getDataName(type, "Amplitudes"))) {
    allocatedTensorArgument<F>(
      getDataName(type, "Amplitudes"), new Tensor<F>(mixer->getNext())
    );
  }
}


template <typename F>
void ClusterSinglesDoublesAlgorithm::amplitudesFromResiduum(
  CTF::Tensor<F> &R, const std::string &indices
) {
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));

  // convert to type F (either complex or double)
  Tensor<F> Fepsi(1, &epsi->lens[0], epsi->sym, *epsi->wrld, "Fepsi");
  // NOTE: just copies if both arguments are real
  toComplexTensor(*epsi, Fepsi);
  Tensor<F> Fepsa(1, &epsa->lens[0], epsa->sym, *epsa->wrld, "Fepsa");
  toComplexTensor(*epsa, Fepsa);

  // create excitation energy tensor
  Tensor<F> D(false, R);
  D.set_name("D");
  for (int p(0); p < D.order/2; ++p) {
    D[indices.c_str()] += Fepsi["i" + p];
    D[indices.c_str()] -= Fepsa["a" + p];
  }

  // TODO:
  // levelshifting can be implemented here

  // use transform to divide T by D
  CTF::Transform<F, F>(
    std::function<void(F, F &)>(
      [](F d, F &t) {
        t = t / d;
      }
    )
  ) (
    D[indices.c_str()], R[indices.c_str()]
  );
}

// instantiate
template
void ClusterSinglesDoublesAlgorithm::amplitudesFromResiduum(
  CTF::Tensor<double> &R, const std::string &indices
);

template
void ClusterSinglesDoublesAlgorithm::amplitudesFromResiduum(
  CTF::Tensor<complex> &R, const std::string &indices
);


template <typename F>
void ClusterSinglesDoublesAlgorithm::dryAmplitudesFromResiduum(
  cc4s::DryTensor<F> &R
) {
  // Build D
  DryTensor<F> D(R, SOURCE_LOCATION);
}

// instantiate
template
void ClusterSinglesDoublesAlgorithm::dryAmplitudesFromResiduum(
  cc4s::DryTensor<double> &R
);

template
void ClusterSinglesDoublesAlgorithm::dryAmplitudesFromResiduum(
  cc4s::DryTensor<complex> &R
);


std::string ClusterSinglesDoublesAlgorithm::getCapitalizedAbbreviation() {
  std::string capitalizedAbbreviation(getAbbreviation());
  std::transform(
    capitalizedAbbreviation.begin(), capitalizedAbbreviation.end(),
    capitalizedAbbreviation.begin(), ::toupper
  );
  return capitalizedAbbreviation;
}


std::string ClusterSinglesDoublesAlgorithm::getDataName(
  const std::string &type, const std::string &data
) {
  std::stringstream dataName;
  dataName << getAbbreviation() << type << data;
  return dataName.str();
}

