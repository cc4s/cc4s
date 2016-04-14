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

/**
 * \brief Calculates the energy of a ClusterDoubles algorithm
 */
void ClusterDoublesAlgorithm::run() {
  // Read the Coulomb Integrals Vabij required for the energy
  Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));

  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  
  // Compute the No,Nv
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // instantiate mixer for the doubles amplitudes, by default use the linear one
  std::string mixerName(getTextArgument("mixer", "LinearMixer"));
  TabijMixer = MixerFactory<double>::create(mixerName, this);
  if (!TabijMixer) {
    std::stringstream stringStream;
    stringStream << "Mixer not implemented: " << mixerName;
    throw new Exception(stringStream.str());
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
  LOG(0, abbreviation) << "Solving Doubles Amplitude Equations" << std::endl;

  // Iteration for determining the DCD amplitudes Tabij
  // and the Dcd energy e
  int64_t maxIterationsCount(
    getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS)
  );
  for (int i(0); i < maxIterationsCount; ++i) {
    LOG(0, abbreviation) << "iteration: " << i+1 << std::endl;
    // call the iterate of the actual algorithm, which is still left open here
    iterate(i);
    energy[""] = 2.0 * TabijMixer->getNext()["abij"] * (*Vabij)["abij"];
    dire = energy.get_val();
    energy[""] = TabijMixer->getNext()["abji"] * (*Vabij)["abij"];
    exce = -1.0 * energy.get_val();
    e = dire + exce;
    LOG(0, abbreviation) << "e=" << e << std::endl;
    LOG(1, abbreviation) << "dir=" << dire << std::endl;
    LOG(1, abbreviation) << "exc=" << exce << std::endl;
  }

  // TODO: use "=" (without space) as soon as most postprocesing
  // is based on TensorWriter
  LOG(1, abbreviation) <<
    abbreviation << " correlation energy = " << e << std::endl;

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
//  DryTensor<> *Vabij(
    getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals");
//  );
  // Read the Particle/Hole Eigenenergies epsi epsa required for the energy
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies")
  );
  
  // Compute the No,Nv
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  std::string abbreviation(getAbbreviation());
  std::transform(
    abbreviation.begin(), abbreviation.end(), abbreviation.begin(),
    ::toupper
  );
  std::stringstream amplitudesName;
  amplitudesName << getAbbreviation() << "DoublesAmplitudes";

  // instantiate mixer for the doubles amplitudes, by default use the linear one
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
    int syms[] = { NS, NS, NS, NS };
    int vvoo[] = { Nv, Nv, No, No };
    DryTensor<> Tabij(4, vvoo, syms);
    allocatedTensorArgument(amplitudesName.str(), new DryTensor<>(Tabij));
  }

  // Allocate the energy e
  DryScalar<> energy();

  LOG(0, abbreviation) << "Solving Doubles Amplitude Equations" << std::endl;

  // Iteration for determining the DCD amplitudes Tabij
  // and the Dcd energy e
//  int64_t maxIterationsCount(
    getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS);
//  );

  // call the dry iterate of the actual algorithm, which is left open here
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


Tensor<> *ClusterDoublesAlgorithm::sliceCoulombIntegrals(int a, int b) {
  Tensor<complex> *GammaGqr(getTensorArgument<complex>("CoulombVertex"));
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(No+Nv);
  int NG(GammaGqr->lens[0]);
  
  // slice the respective parts from the Coulomb vertex
  int leftGammaStart[] = { 0, No+a, No };
  int leftGammaEnd[] = { NG, std::min(No+a+No, Np), Np };
  int rightGammaStart[] = { 0, No+b, No };
  int rightGammaEnd[] = { NG, std::min(No+b+No, Np), Np };
  Tensor<complex> leftGamma(GammaGqr->slice(leftGammaStart, leftGammaEnd));
  Tensor<complex> rightGamma(GammaGqr->slice(rightGammaStart, rightGammaEnd));
  // split into real and imaginary parts
  Tensor<> realLeftGamma(3, leftGamma.lens, leftGamma.sym, *GammaGqr->wrld);
  Tensor<> imagLeftGamma(3, leftGamma.lens, leftGamma.sym, *GammaGqr->wrld);
  fromComplexTensor(leftGamma, realLeftGamma, imagLeftGamma);
  Tensor<> realRightGamma(3, rightGamma.lens, rightGamma.sym, *GammaGqr->wrld);
  Tensor<> imagRightGamma(3, rightGamma.lens, rightGamma.sym, *GammaGqr->wrld);
  fromComplexTensor(rightGamma, realRightGamma, imagRightGamma);

  // allocate sliced Coulomb integrals
  int lens[] = {
    leftGamma.lens[1], rightGamma.lens[1], leftGamma.lens[2], rightGamma.lens[2]
  };
  int syms[] = { NS, NS, NS, NS };
  Tensor<> *Vxycd(new Tensor<>(4, lens, syms, *GammaGqr->wrld));

  // contract left and right slices of the Coulomb vertices
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
    // add the same slice at (b,a,j,i):
    dstStart[0] = b; dstStart[1] = a;
    dstEnd[0] = b+Ny; dstEnd[1] = a+Nx;
    srcEnd[0] = Ny; srcEnd[1] = Nx;
    // swap xy and ij simultaneously
    Tensor<> Ryxji(4, srcEnd, Rxyij.sym, *Rxyij.wrld);
    Ryxji["yxji"] = Rxyij["xyij"];
    // and add it
    Rabij.slice(dstStart,dstEnd,1.0, Ryxji,srcStart,srcEnd,1.0);
  }
}

