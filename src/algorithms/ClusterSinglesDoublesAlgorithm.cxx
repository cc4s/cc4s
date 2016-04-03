#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <Options.hpp>

using namespace CTF;
using namespace cc4s;

ClusterSinglesDoublesAlgorithm::ClusterSinglesDoublesAlgorithm(
  std::vector<Argument> const &argumentList
): ClusterDoublesAlgorithm(argumentList) {
}

ClusterSinglesDoublesAlgorithm::~ClusterSinglesDoublesAlgorithm() {
}

/**
 * \brief Calculates the energy of a ClusterSinglesDoubles algorithm
 */
void ClusterSinglesDoublesAlgorithm::run() {
  // Read the Coulomb Integrals Vabij required for the energy
  Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  
  // Compute the No,Nv
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // TODO: factor out code common with ClusterDoublesAlgorithm
  // instantiate mixer for the doubles amplitudes, by default use the linear one
  std::string mixerName(getTextArgument("mixer", "LinearMixer"));
  TabijMixer = MixerFactory<double>::create(mixerName, this);
  if (!TabijMixer) {
    std::stringstream stringStream;
    stringStream << "Mixer not implemented: " << mixerName;
    throw new Exception(stringStream.str());
  }
  // create another mixer for the singles, this time it must exists
  TaiMixer = MixerFactory<double>::create(mixerName, this);
  {
    // Allocate the singles amplitudes and append it to the mixer
    int syms[] = { NS, NS };
    int vo[] = { Nv, No };
    Tensor<> Tai(2, vo, syms, *epsi->wrld, "Tai");
    TaiMixer->append(Tai);
    // the amplitudes will from now on be managed by the mixer
  }
  {
    // Allocate the doubles amplitudes and append it to the mixer
    int syms[] = { NS, NS, NS, NS };
    int vvoo[] = { Nv, Nv, No, No };
    Tensor<> Tabij(4, vvoo, syms, *epsi->wrld, "Tabij");
    TabijMixer->append(Tabij);
    // the amplitudes will from now on be managed by the mixer
  }

  // Allocate the energy e
  Scalar<> energy(*epsi->wrld);
  double e(0), dire, exce;

  std::string abbreviation(getAbbreviation());
  std::transform(
    abbreviation.begin(), abbreviation.end(), abbreviation.begin(),
    ::toupper
  );
  LOG(0, abbreviation) << "Solving Singles and Doubles Amplitude Equations" << std::endl;

  // Iteration for determining the DCD amplitudes Tabij
  // and the Dcd energy e
  int64_t maxIterationsCount(
    getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS)
  );
  for (int i(0); i < maxIterationsCount; ++i) {
    LOG(0, abbreviation) << "iteration: " << i+1 << std::endl;
    // call the iterate of the actual algorithm, which is still left open here
    iterate(i);
    Tensor<> *Tai(&TaiMixer->getNext());
    Tensor<> *Tabij(&TabijMixer->getNext());
    // Singles direct term
    energy[""]  = 2.0 * (*Vabij)["abij"] * (*Tai)["ai"] * (*Tai)["bj"];
    // Doubles direct term
    energy[""] += 2.0 * (*Tabij)["abij"] * (*Vabij)["abij"];
    // Compute direct energy
    dire = energy.get_val();
    // Singles exchange term
    energy[""]  = (*Vabij)["baij"] * (*Tai)["ai"] * (*Tai)["bj"];
    // Doubles exchange term
    energy[""] += (*Tabij)["abji"] * (*Vabij)["abij"];
    // Compute exchange energy
    exce = -1.0 * energy.get_val();
    // Compute total energy
    e = dire + exce;
    LOG(0, abbreviation) << "e=" << e << std::endl;
    LOG(1, abbreviation) << "dir=" << dire << std::endl;
    LOG(1, abbreviation) << "exc=" << exce << std::endl;
  }

  // TODO: use "=" (without space) as soon as most postprocesing
  // is based on TensorWriter
  LOG(1, abbreviation) <<
    abbreviation << " correlation energy = " << e << std::endl;

  std::stringstream doublesAmplitudesName;
  doublesAmplitudesName << getAbbreviation() << "DoublesAmplitudes";
  allocatedTensorArgument(
    doublesAmplitudesName.str(), new Tensor<>(TabijMixer->getNext())
  );

  std::stringstream singlesAmplitudesName;
  singlesAmplitudesName << getAbbreviation() << "SinglesAmplitudes";
  allocatedTensorArgument(
    singlesAmplitudesName.str(), new Tensor<>(TaiMixer->getNext())
  );

  std::stringstream energyName;
  energyName << getAbbreviation() << "Energy";
  setRealArgument(energyName.str(), e);
}

void ClusterSinglesDoublesAlgorithm::singlesAmplitudesFromResiduum(
  CTF::Tensor<> &Rai
) {
  // Build Dai
  Tensor<> Dai(false, Rai);
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  Dai["ai"]  = (*epsi)["i"];
  Dai["ai"] -= (*epsa)["a"];

  // TODO:
  // levelshifting can be implemented here...

  // Divide Rai/Dai to get Tai
  Bivar_Function<> fDivide(&divide<double>);
  Rai.contract(1.0, Rai,"ai", Dai,"ai", 0.0,"ai", fDivide);
}

Tensor<> *ClusterSinglesDoublesAlgorithm::sliceCoupledCoulombIntegrals(
  int a, int b
) {
  // get the sliced Coulomb integrals Vxycd
  Tensor<> *Xxycd(sliceCoulombIntegrals(a, b));

  // couple to singles
  Tensor<> *Vabci(getTensorArgument("PPPHCoulombIntegrals"));
  Tensor<> *Tai(&TaiMixer->getNext());
  int Nx(Xxycd->lens[0]);
  int Ny(Xxycd->lens[1]);
  int No(Tai->lens[1]);
  int Nv(Tai->lens[0]);
  // calculate Xxycd["xycd"] -= (*Vabci)["cdxk"] * (*Tai)["yk"];
  int VStart[] = { 0, 0, a, 0 }; int VEnd[] = { Nv, Nv, a+Nx, No };
  int TStart[] = { b, 0 }; int TEnd[] = { b+Ny, No };
  (*Xxycd)["xycd"] -=
    Vabci->slice(VStart,VEnd)["cdxk"] * Tai->slice(TStart,TEnd)["yk"];
  // calculate Xxycd["xycd"] -= (*Vabyi)["dcyk"] * (*Txi)["xk"];
  VStart[2] = b; VEnd[2] = b+Ny;
  TStart[0] = a; TEnd[0] = a+Nx;
  (*Xxycd)["xycd"] -=
    Vabci->slice(VStart,VEnd)["dcyk"] * Tai->slice(TStart,TEnd)["xk"];

  return Xxycd;
}

