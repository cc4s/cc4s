#include <algorithms/ClusterSinglesDoublesAlgorithm.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <util/DryTensor.hpp>
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
  if (TaiMixer) delete TaiMixer;
}

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
    // The amplitudes will from now on be managed by the mixer
  }

  // Allocate the energy e
  Scalar<> energy(*epsi->wrld);
  energy.set_name("energy");
  double e(0), dire, exce;

  std::string abbreviation(getAbbreviation());
  std::transform(abbreviation.begin(), abbreviation.end(), 
		 abbreviation.begin(), ::toupper);

  // Iteration for determining the amplitudes Tai and Tabij
  // and the energy e
  int maxIterationsCount(getIntegerArgument("maxIterations", 
						DEFAULT_MAX_ITERATIONS));

  for (int i(0); i < maxIterationsCount; ++i) {
    LOG(0, abbreviation) << "iteration: " << i+1 << std::endl;
    // Call the iterate of the actual algorithm, which is still left open here
    iterate(i);
    Tensor<> *Tai(&TaiMixer->getNext());
    Tai->set_name("Tai");
    Tensor<> *Tabij(&TabijMixer->getNext());
    Tabij->set_name("Tabij");
    // Singles direct term
    energy[""]  = 2.0 * (*Tai)["bj"] * (*Vabij)["abij"] * (*Tai)["ai"];
    // Doubles direct term
    energy[""] += 2.0 * (*Tabij)["abij"] * (*Vabij)["abij"];
    // Compute direct energy
    dire = energy.get_val();
    // Singles exchange term
    energy[""]  =  (*Tai)["bj"] * (*Vabij)["baij"] * (*Tai)["ai"];
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

void ClusterSinglesDoublesAlgorithm::dryRun() {
  // Read the Coulomb Integrals Vabij required for the energy
  getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals");

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(getTensorArgument<double, 
		    DryTensor<double>>("HoleEigenEnergies"));
  DryTensor<> *epsa(getTensorArgument<double, 
		    DryTensor<double>>("ParticleEigenEnergies"));

  std::string abbreviation(getAbbreviation());
  std::transform(abbreviation.begin(), abbreviation.end(), 
		 abbreviation.begin(), ::toupper);

  // Instantiate mixer for the doubles amplitudes, by default use the linear one
  std::string mixerName(getTextArgument("mixer", "LinearMixer"));
  TabijMixer = MixerFactory<double>::create(mixerName, this);
  if (!TabijMixer) {
    std::stringstream stringStream;
    stringStream << "Mixer not implemented: " << mixerName;
    throw new Exception(stringStream.str());
  }
  // TODO: implement DryTensor in mixers
  if (mixerName != "LinearMixer") {
    LOG(0, abbreviation)
      << "Warning: dry run not implemented for " << mixerName
      << ", assuming the same memory usage." << std::endl;
  }

  {
    // Allocate the doubles amplitudes and append it to the mixer
    int No(epsi->lens[0]);
    int Nv(epsa->lens[0]);
    int syms[] = { NS, NS, NS, NS };
    int vvoo[] = { Nv, Nv, No, No };
    int vo[] = { Nv, No };
    std::stringstream doublesAmplitudesName;
    doublesAmplitudesName << getAbbreviation() << "DoublesAmplitudes";
    DryTensor<> Tabij(4, vvoo, syms);
    allocatedTensorArgument(doublesAmplitudesName.str(), 
			    new DryTensor<>(Tabij));

    std::stringstream singlesAmplitudesName;
    singlesAmplitudesName << getAbbreviation() << "SinglesAmplitudes";
    DryTensor<> Tai(2, vo, syms);
    allocatedTensorArgument(singlesAmplitudesName.str(), 
			    new DryTensor<>(Tai));
  }

  // Allocate the energy e
  DryScalar<> energy();

  getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS);

  // Call the dry iterate of the actual algorithm, which is left open here
  dryIterate();

  std::stringstream energyName;
  energyName << getAbbreviation() << "Energy";
  setRealArgument(energyName.str(), 0.0);
}

void ClusterSinglesDoublesAlgorithm::singlesAmplitudesFromResiduum(
  CTF::Tensor<> &Rai
) {
  // Build Dai
  Tensor<> Dai(false, Rai);
  Dai.set_name("Dai");
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
  int a, int b, int sliceRank) 
{
  // Get the sliced Coulomb integrals Vxycd 
  Tensor<> *Xxycd(sliceCoulombIntegrals(a, b, sliceRank));
  Xxycd->set_name("Xxycd");

  // Couple to singles
  Tensor<> *Vabci(getTensorArgument("PPPHCoulombIntegrals"));
  Tensor<> *Tai(&TaiMixer->getNext());
  int Nx(Xxycd->lens[0]);
  int Ny(Xxycd->lens[1]);
  int No(Tai->lens[1]);
  int Nv(Tai->lens[0]);

  // Calculate Xxycd["xycd"] -= (*Vabci)["cdxk"] * (*Tai)["yk"];
  int VStart[] = { 0, 0, a, 0 }; int VEnd[] = { Nv, Nv, a+Nx, No };
  int TStart[] = { b, 0 }; int TEnd[] = { b+Ny, No };
  (*Xxycd)["xycd"] -=
    Vabci->slice(VStart,VEnd)["cdxk"] * Tai->slice(TStart,TEnd)["yk"];

  // Calculate Xxycd["xycd"] -= (*Vabyi)["dcyk"] * (*Txi)["xk"];
  VStart[2] = b; VEnd[2] = b+Ny;
  TStart[0] = a; TEnd[0] = a+Nx;
  (*Xxycd)["xycd"] -=
    Vabci->slice(VStart,VEnd)["dcyk"] * Tai->slice(TStart,TEnd)["xk"];

  return Xxycd;
}

