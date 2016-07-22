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
  
  // Compute the No,Nv
  int No(Vabij->lens[2]);
  int Nv(Vabij->lens[0]);

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
    Tensor<> Tai(2, vo, syms, *Vabij->wrld, "Tai");
    TaiMixer->append(Tai);
    // the amplitudes will from now on be managed by the mixer
  }
  {
    // Allocate the doubles amplitudes and append it to the mixer
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

    // Intermediate tensor Xabij=T2+T1*T1
    Tensor<> Xabij(Tabij);
    Xabij.set_name("Xabij");
    Xabij["abij"] += (*Tai)["ai"] * (*Tai)["bj"];

    // Direct term
    energy[""] = 2.0 * Xabij["abij"] * (*Vabij)["abij"];
    // Compute direct energy
    dire = energy.get_val();
    // Doubles exchange term
    energy[""] = Xabij["abij"] * (*Vabij)["baij"];
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

void ClusterSinglesDoublesAlgorithm::drySinglesAmplitudesFromResiduum(
  cc4s::DryTensor<> &Rai
) {
  // Build Dai
  DryTensor<> Dai(Rai);
}

Tensor<> *ClusterSinglesDoublesAlgorithm::sliceCoupledCoulombIntegrals(int a, 
								       int b, 
								       int sliceRank)
{
  // Read the amplitudes Tai
  Tensor<> *Tai(&TaiMixer->getNext());
  Tai->set_name("Tai");

  // Read the Coulomb vertex GammaGpq
  Tensor<complex> *GammaGpq( getTensorArgument<complex>("CoulombVertex"));
  GammaGpq->set_name("GammaGpq");

  // Compute No,Nv,NG,Np
  int No(Tai->lens[1]);
  int Nv(Tai->lens[0]);
  int NG(GammaGpq->lens[0]);
  int Np = No + Nv;

  // Allocate and compute GammaGab,GammaGai from GammaGpq
  int GaiStart[] = {0 ,No, 0};
  int GaiEnd[]   = {NG,Np,No};
  int GabStart[] = {0 ,No,No};
  int GabEnd[]   = {NG,Np,Np};
  Tensor<complex> GammaGai(GammaGpq->slice(GaiStart,GaiEnd));
  Tensor<complex> GammaGab(GammaGpq->slice(GabStart,GabEnd));

  // Split GammaGab,GammaGai into real and imaginary parts
  Tensor<> realGammaGai(3, GammaGai.lens, GammaGai.sym, 
			*GammaGai.wrld, "RealGammaGai");
  Tensor<> imagGammaGai(3, GammaGai.lens, GammaGai.sym, 
			*GammaGai.wrld, "ImagGammaGai");
  fromComplexTensor(GammaGai, realGammaGai, imagGammaGai);

  Tensor<> realGammaGab(3, GammaGab.lens, GammaGab.sym, 
			*GammaGab.wrld, "RealGammaGab");
  Tensor<> imagGammaGab(3, GammaGab.lens, GammaGab.sym, 
			*GammaGab.wrld, "ImagGammaGab");
  fromComplexTensor(GammaGab, realGammaGab, imagGammaGab);

  // Construct dressed Coulomb vertex GammaGab
  realGammaGab["Gab"] += (-1.0) * realGammaGai["Gbk"] * (*Tai)["ak"];
  imagGammaGab["Gab"] += (-1.0) * imagGammaGai["Gbk"] * (*Tai)["ak"];
  toComplexTensor(realGammaGab, imagGammaGab, GammaGab);
  
  // Slice the respective parts from the dressed Coulomb vertex GammaGab
  int leftGammaStart[] = { 0, a, 0 };
  int leftGammaEnd[] = { NG, std::min(a+sliceRank, Nv), Nv };
  int rightGammaStart[] = { 0, b, 0 };
  int rightGammaEnd[] = { NG, std::min(b+sliceRank, Nv), Nv };
  
  Tensor<complex> leftGamma(GammaGab.slice(leftGammaStart, leftGammaEnd));
  Tensor<complex> rightGamma(GammaGab.slice(rightGammaStart, rightGammaEnd));

  // Split into real and imaginary parts
  Tensor<> realLeftGamma(3, leftGamma.lens, leftGamma.sym, *GammaGpq->wrld, "realLeftGamma");
  Tensor<> imagLeftGamma(3, leftGamma.lens, leftGamma.sym, *GammaGpq->wrld, "imagLeftGamma");
  fromComplexTensor(leftGamma, realLeftGamma, imagLeftGamma);
  Tensor<> realRightGamma(3, rightGamma.lens, rightGamma.sym, *GammaGpq->wrld, "realRightGamma");
  Tensor<> imagRightGamma(3, rightGamma.lens, rightGamma.sym, *GammaGpq->wrld, "imagRightGamma");
  fromComplexTensor(rightGamma, realRightGamma, imagRightGamma);
  
  // Allocate sliced Coulomb integrals
  int lens[] = {
    leftGamma.lens[1], rightGamma.lens[1], leftGamma.lens[2], rightGamma.lens[2]
  };
  int syms[] = { NS, NS, NS, NS };
  Tensor<> *Vxycd(new Tensor<>(4, lens, syms, *GammaGpq->wrld, "Vxycd"));

  // Contract left and right slices of the dressed Coulomb vertices
  (*Vxycd)["xycd"]  = realLeftGamma["Gxc"] * realRightGamma["Gyd"];
  (*Vxycd)["xycd"] += imagLeftGamma["Gxc"] * imagRightGamma["Gyd"];
  return Vxycd;
}

DryTensor<> *ClusterSinglesDoublesAlgorithm::drySliceCoupledCoulombIntegrals(int sliceRank)
{
  // Read the Coulomb vertex GammaGpq
  DryTensor<complex> *GammaGpq(getTensorArgument<complex, 
			       DryTensor<complex>>("CoulombVertex"));
  
  // Read the Particle/Hole Eigenenergies
  DryTensor<> *epsi(getTensorArgument
		    <double, DryTensor<double>>("HoleEigenEnergies"));
  DryTensor<> *epsa(getTensorArgument
		    <double, DryTensor<double>>("ParticleEigenEnergies"));

  // Compute No,Nv,NG,Np
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int NG(GammaGpq->lens[0]);
  int syms[] = { NS, NS, NS, NS };

  // Allocate and compute GammaGab,GammaGai from GammaGpq
  int GaiLens[]   = {NG,Nv,No};
  int GabLens[]   = {NG,Nv,Nv};
  int GijLens[]   = {NG,No,No};

  DryTensor<complex> GammaGai(3, GaiLens, syms);
  DryTensor<complex> GammaGab(3, GabLens, syms);
  DryTensor<complex> GammaGij(3, GijLens, syms);

  // Split GammaGab,GammaGai into real and imaginary parts
  DryTensor<> realGammaGai(3, GaiLens, syms);
  DryTensor<> imagGammaGai(3, GaiLens, syms);

  DryTensor<> realGammaGab(3, GabLens, syms);
  DryTensor<> imagGammaGab(3, GabLens, syms);

  DryTensor<> realGammaGij(3, GijLens, syms);
  DryTensor<> imagGammaGij(3, GijLens, syms);
  
  // Slice the respective parts from the Coulomb vertex
  int leftGammaLens[]  = { NG, sliceRank, Nv };
  int rightGammaLens[] = { NG, sliceRank, Nv };
  DryTensor<complex> leftGamma (3, leftGammaLens , syms);
  DryTensor<complex> rightGamma(3, rightGammaLens, syms);

  // Split into real and imaginary parts
  DryTensor<> realLeftGamma(3, leftGammaLens, syms);
  DryTensor<> imagLeftGamma(3, leftGammaLens, syms);

  DryTensor<> realRightGamma(3, rightGammaLens, syms);
  DryTensor<> imagRightGamma(3, rightGammaLens, syms);

  // Allocate sliced Coulomb integrals
  int lens[] = {leftGamma.lens[1], rightGamma.lens[1], 
		leftGamma.lens[2], rightGamma.lens[2]};
  DryTensor<> *Vxycd(new DryTensor<>(4, lens, syms));

  return Vxycd;
}

