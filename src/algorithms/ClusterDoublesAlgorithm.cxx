#include <algorithms/ClusterDoublesAlgorithm.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <util/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <Options.hpp>

using namespace CTF;
using namespace cc4s;

ClusterDoublesAlgorithm::ClusterDoublesAlgorithm(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ClusterDoublesAlgorithm::~ClusterDoublesAlgorithm() {
  if (TabijMixer) delete TabijMixer;
}

void ClusterDoublesAlgorithm::run() {
  // Read the Coulomb Integrals Vabij required for the energy
  Tensor<> *Vabij(getTensorArgument<>("PPHHCoulombIntegrals"));

  // Instantiate mixer for the doubles amplitudes, by default use the linear one
  std::string mixerName(getTextArgument("mixer", "LinearMixer"));
  TabijMixer = MixerFactory<double>::create(mixerName, this);
  if (!TabijMixer) {
    std::stringstream stringStream;
    stringStream << "Mixer not implemented: " << mixerName;
    throw new Exception(stringStream.str());
  }

  {
    // Allocate the doubles amplitudes and append it to the mixer
    int No(Vabij->lens[2]);
    int Nv(Vabij->lens[0]);
    int syms[] = { NS, NS, NS, NS };
    int vvoo[] = { Nv, Nv, No, No };
    Tensor<> Tabij(4, vvoo, syms, *Vabij->wrld, "Tabij");
    TabijMixer->append(Tabij);
    // The amplitudes will from now on be managed by the mixer
  }

  // Allocate the energy e
  Scalar<> energy(*Vabij->wrld);
  energy.set_name("energy");
  double e(0), dire, exce;

  std::string abbreviation(getAbbreviation());
  std::transform(abbreviation.begin(), abbreviation.end(), 
		 abbreviation.begin(), ::toupper);

  // Iteration for determining the doubles amplitudes Tabij
  // and the energy e
  int maxIterationsCount(getIntegerArgument("maxIterations", 
					    DEFAULT_MAX_ITERATIONS));
  for (int i(0); i < maxIterationsCount; ++i) {
    LOG(0, abbreviation) << "iteration: " << i+1 << std::endl;
    // call the iterate of the actual algorithm, which is still left open here
    iterate(i);
    Tensor<> *Tabij(&TabijMixer->getNext());
    Tabij->set_name("Tabij");
    // Direct term
    energy[""] = 2.0 * (*Tabij)["abij"] * (*Vabij)["abij"];
    dire = energy.get_val();
    // Exchange term
    energy[""] = (*Tabij)["abji"] * (*Vabij)["abij"];
    exce = -1.0 * energy.get_val();
    // Total energy
    e = dire + exce;
    LOG(0, abbreviation) << "e=" << e << std::endl;
    LOG(1, abbreviation) << "dir=" << dire << std::endl;
    LOG(1, abbreviation) << "exc=" << exce << std::endl;
  }

  std::stringstream amplitudesName;
  amplitudesName << getAbbreviation() << "DoublesAmplitudes";
  allocatedTensorArgument(
    amplitudesName.str(), new Tensor<>(TabijMixer->getNext())
  );

  std::stringstream energyName;
  energyName << getAbbreviation() << "Energy";
  setRealArgument(energyName.str(), e);
}

void ClusterDoublesAlgorithm::dryRun() {
  // Read the Coulomb Integrals Vabij required for the energy
  getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals");

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(getTensorArgument<double, 
		    DryTensor<double>>("HoleEigenEnergies"));
  DryTensor<> *epsa(getTensorArgument<double, 
		    DryTensor<double>>("ParticleEigenEnergies"));

  std::string abbreviation(getAbbreviation());
  std::transform(
    abbreviation.begin(), abbreviation.end(), abbreviation.begin(),
    ::toupper
  );
  std::stringstream amplitudesName;
  amplitudesName << getAbbreviation() << "DoublesAmplitudes";

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
    DryTensor<> Tabij(4, vvoo, syms);
    allocatedTensorArgument(amplitudesName.str(), new DryTensor<>(Tabij));
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

void ClusterDoublesAlgorithm::dryIterate() {
  LOG(0, "CluserDoubles") << "Dry run for iteration not given for "
    << getAbbreviation() << std::endl;
}

void ClusterDoublesAlgorithm::doublesAmplitudesFromResiduum(
  CTF::Tensor<> &Rabij
) {
  // Build Dabij
  Tensor<> Dabij(false, Rabij);
  Dabij.set_name("Dabij");
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  Dabij["abij"]  = (*epsi)["i"];
  Dabij["abij"] += (*epsi)["j"];
  Dabij["abij"] -= (*epsa)["a"];
  Dabij["abij"] -= (*epsa)["b"];

  // TODO:
  // levelshifting can be implemented here

  // Divide Rabij/Dabij to get Tabij
  Bivar_Function<> fDivide(&divide<double>);
  Rabij.contract(1.0, Rabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);
}


Tensor<> *ClusterDoublesAlgorithm::sliceCoulombIntegrals(int a, int b, int sliceRank) {
  Tensor<complex> *GammaGqr(getTensorArgument<complex>("CoulombVertex"));
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(No+Nv);
  int NG(GammaGqr->lens[0]);
  
  // Slice the respective parts from the Coulomb vertex
  int leftGammaStart[] = { 0, No+a, No };
  int leftGammaEnd[] = { NG, std::min(No+a+sliceRank, Np), Np };
  int rightGammaStart[] = { 0, No+b, No };
  int rightGammaEnd[] = { NG, std::min(No+b+sliceRank, Np), Np };
  Tensor<complex> leftGamma(GammaGqr->slice(leftGammaStart, leftGammaEnd));
  Tensor<complex> rightGamma(GammaGqr->slice(rightGammaStart, rightGammaEnd));
  // Split into real and imaginary parts
  Tensor<> realLeftGamma(3, leftGamma.lens, leftGamma.sym, *GammaGqr->wrld, "realLeftGamma");
  Tensor<> imagLeftGamma(3, leftGamma.lens, leftGamma.sym, *GammaGqr->wrld, "imagLeftGamma");
  fromComplexTensor(leftGamma, realLeftGamma, imagLeftGamma);
  Tensor<> realRightGamma(3, rightGamma.lens, rightGamma.sym, *GammaGqr->wrld, "realRightGamma");
  Tensor<> imagRightGamma(3, rightGamma.lens, rightGamma.sym, *GammaGqr->wrld, "imagRightGamma");
  fromComplexTensor(rightGamma, realRightGamma, imagRightGamma);

  // Allocate sliced Coulomb integrals
  int lens[] = {
    leftGamma.lens[1], rightGamma.lens[1], leftGamma.lens[2], rightGamma.lens[2]
  };
  int syms[] = { NS, NS, NS, NS };
  Tensor<> *Vxycd(new Tensor<>(4, lens, syms, *GammaGqr->wrld, "Vxycd"));

  // Contract left and right slices of the Coulomb vertices
  (*Vxycd)["xycd"] =  realLeftGamma["Gxc"] * realRightGamma["Gyd"];
  (*Vxycd)["xycd"] += imagLeftGamma["Gxc"] * imagRightGamma["Gyd"];
  return Vxycd;
}

void ClusterDoublesAlgorithm::sliceIntoResiduum(
  Tensor<> &Rxyij, int a, int b, Tensor<> &Rabij
) {
  int Nx(Rxyij.lens[0]);
  int Ny(Rxyij.lens[1]);
  int No(Rxyij.lens[2]);
  int dstStart[] = { a, b, 0, 0 };
  int dstEnd[] = { a+Nx, b+Ny, No, No };
  int srcStart[] = { 0, 0, 0, 0 };
  int srcEnd[] = { Nx, Ny, No, No };
  // R["abij"] += R["xyij"] at current x,y
  Rabij.slice(dstStart,dstEnd,1.0, Rxyij,srcStart,srcEnd,1.0);
  if (a>b) {
    // Add the same slice at (b,a,j,i):
    dstStart[0] = b; dstStart[1] = a;
    dstEnd[0] = b+Ny; dstEnd[1] = a+Nx;
    srcEnd[0] = Ny; srcEnd[1] = Nx;
    // Swap xy and ij simultaneously
    Tensor<> Ryxji(4, srcEnd, Rxyij.sym, *Rxyij.wrld, "Ryxji");
    Ryxji["yxji"] = Rxyij["xyij"];
    // Add Ryxij to Rabij
    Rabij.slice(dstStart,dstEnd,1.0, Ryxji,srcStart,srcEnd,1.0);
  }
}

