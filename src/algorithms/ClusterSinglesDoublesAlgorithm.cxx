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
    if (isArgumentGiven("StartingDoublesAmplitudes")) {
      Tensor<> Tabij(getTensorArgument("StartingDoublesAmplitudes"));
      TabijMixer->append(Tabij);
    }
    else {
      Tensor<> Tabij(4, vvoo, syms, *Vabij->wrld, "Tabij");
      TabijMixer->append(Tabij);
    }
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
  if (isArgumentGiven(doublesAmplitudesName.str())) {
    allocatedTensorArgument(
      doublesAmplitudesName.str(), new Tensor<>(TabijMixer->getNext())
    );
  }

  std::stringstream singlesAmplitudesName;
  singlesAmplitudesName << getAbbreviation() << "SinglesAmplitudes";
  if (isArgumentGiven(singlesAmplitudesName.str())) {
    allocatedTensorArgument(
      singlesAmplitudesName.str(), new Tensor<>(TaiMixer->getNext())
    );
  }

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
    DryTensor<> Tabij(4, vvoo, syms, SOURCE_LOCATION);
    allocatedTensorArgument(doublesAmplitudesName.str(), 
                            new DryTensor<>(Tabij, SOURCE_LOCATION));

    std::stringstream singlesAmplitudesName;
    singlesAmplitudesName << getAbbreviation() << "SinglesAmplitudes";
    DryTensor<> Tai(2, vo, syms, SOURCE_LOCATION);
    allocatedTensorArgument(singlesAmplitudesName.str(), 
                            new DryTensor<>(Tai, SOURCE_LOCATION));
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
  DryTensor<> Dai(Rai, SOURCE_LOCATION);
}

Tensor<> *ClusterSinglesDoublesAlgorithm::sliceCoupledCoulombIntegrals(int a, 
                                                                       int b, 
                                                                       int integralsSliceSize)
{
  // Read the amplitudes Tai
  Tensor<> *Tai(&TaiMixer->getNext());
  Tai->set_name("Tai");

  // Read the Coulomb vertex GammaGqr
  Tensor<complex> *GammaGqr( getTensorArgument<complex>("CoulombVertex"));
  GammaGqr->set_name("GammaGqr");

  // Compute No,Nv,NG,Np
  int No(Tai->lens[1]);
  int Nv(Tai->lens[0]);
  int NG(GammaGqr->lens[0]);
  int Np = No + Nv;

  // Allocate and compute GammaGab,GammaGai from GammaGqr
  int GaiStart[] = {0 ,No, 0};
  int GaiEnd[]   = {NG,Np,No};
  int GabStart[] = {0 ,No,No};
  int GabEnd[]   = {NG,Np,Np};
  Tensor<complex> GammaGai(GammaGqr->slice(GaiStart,GaiEnd));
  Tensor<complex> GammaGab(GammaGqr->slice(GabStart,GabEnd));

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
  int leftGammaEnd[] = { NG, std::min(a+integralsSliceSize, Nv), Nv };
  int rightGammaStart[] = { 0, b, 0 };
  int rightGammaEnd[] = { NG, std::min(b+integralsSliceSize, Nv), Nv };
  
  Tensor<complex> leftGamma(GammaGab.slice(leftGammaStart, leftGammaEnd));
  Tensor<complex> rightGamma(GammaGab.slice(rightGammaStart, rightGammaEnd));

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

  // Contract left and right slices of the dressed Coulomb vertices
  (*Vxycd)["xycd"]  = realLeftGamma["Gxc"] * realRightGamma["Gyd"];
  (*Vxycd)["xycd"] += imagLeftGamma["Gxc"] * imagRightGamma["Gyd"];
  return Vxycd;
}

DryTensor<> *ClusterSinglesDoublesAlgorithm::drySliceCoupledCoulombIntegrals(int integralsSliceSize)
{
  // Read the Coulomb vertex GammaGqr
  DryTensor<complex> *GammaGqr(getTensorArgument<complex, 
                               DryTensor<complex>>("CoulombVertex"));
  
  // Read the Particle/Hole Eigenenergies
  DryTensor<> *epsi(getTensorArgument
                    <double, DryTensor<double>>("HoleEigenEnergies"));
  DryTensor<> *epsa(getTensorArgument
                    <double, DryTensor<double>>("ParticleEigenEnergies"));

  // Compute No,Nv,NG,Np
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int NG(GammaGqr->lens[0]);
  int syms[] = { NS, NS, NS, NS };

  // Allocate and compute GammaGab,GammaGai from GammaGqr
  int GaiLens[]   = {NG,Nv,No};
  int GabLens[]   = {NG,Nv,Nv};
  int GijLens[]   = {NG,No,No};

  DryTensor<complex> GammaGai(3, GaiLens, syms, SOURCE_LOCATION);
  DryTensor<complex> GammaGab(3, GabLens, syms, SOURCE_LOCATION);
  DryTensor<complex> GammaGij(3, GijLens, syms, SOURCE_LOCATION);

  // Split GammaGab,GammaGai into real and imaginary parts
  DryTensor<> realGammaGai(3, GaiLens, syms, SOURCE_LOCATION);
  DryTensor<> imagGammaGai(3, GaiLens, syms, SOURCE_LOCATION);

  DryTensor<> realGammaGab(3, GabLens, syms, SOURCE_LOCATION);
  DryTensor<> imagGammaGab(3, GabLens, syms, SOURCE_LOCATION);

  DryTensor<> realGammaGij(3, GijLens, syms, SOURCE_LOCATION);
  DryTensor<> imagGammaGij(3, GijLens, syms, SOURCE_LOCATION);
  
  // Slice the respective parts from the Coulomb vertex
  int leftGammaLens[]  = { NG, integralsSliceSize, Nv };
  int rightGammaLens[] = { NG, integralsSliceSize, Nv };
  DryTensor<complex> leftGamma (3, leftGammaLens , syms, SOURCE_LOCATION);
  DryTensor<complex> rightGamma(3, rightGammaLens, syms, SOURCE_LOCATION);

  // Split into real and imaginary parts
  DryTensor<> realLeftGamma(3, leftGammaLens, syms, SOURCE_LOCATION);
  DryTensor<> imagLeftGamma(3, leftGammaLens, syms, SOURCE_LOCATION);

  DryTensor<> realRightGamma(3, rightGammaLens, syms, SOURCE_LOCATION);
  DryTensor<> imagRightGamma(3, rightGammaLens, syms, SOURCE_LOCATION);

  // Allocate sliced Coulomb integrals
  int lens[] = {leftGamma.lens[1], rightGamma.lens[1], 
                leftGamma.lens[2], rightGamma.lens[2]};
  DryTensor<> *Vxycd(new DryTensor<>(4, lens, syms, SOURCE_LOCATION));

  return Vxycd;
}

Tensor<> *ClusterSinglesDoublesAlgorithm::sliceAmplitudesFromCoupledCoulombFactors(int a,
                                                                                   int b,
                                                                                   int factorsSliceSize)
{
  Tensor<complex> *PirR(getTensorArgument<complex>("FactorOrbitals"));
  PirR->set_name("PirR");
  Tensor<complex> *LambdaGR(getTensorArgument<complex>("CoulombFactors"));
  LambdaGR->set_name("LambdaGR");

  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));

  // Read the doubles amplitudes Tabij
  Tensor<> *Tabij(&TabijMixer->getNext());
  Tabij->set_name("Tabij");
  Tensor<> *Tai(&TaiMixer->getNext());
  Tai->set_name("Tai");

  // Intermediate tensor Iabij=T2+T1*T1
  Tensor<> Iabij(Tabij);
  Iabij.set_name("Iabij");
  Iabij["abij"] += (*Tai)["ai"] * (*Tai)["bj"];

  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(PirR->lens[0]);
  int NR(PirR->lens[1]);
  int NG(LambdaGR->lens[0]);
  int Rx(std::min(factorsSliceSize, NR-a));
  int Ry(std::min(factorsSliceSize, NR-b));
  int Rvoo[] = { Rx, Nv, No, No };
  int RRoo[] = { Rx, Ry, No, No };
  int RR[] = { Rx, Ry };
  int syms[] = { NS, NS, NS, NS };

  Tensor<complex> VRS(2, RR, syms, *PirR->wrld, "VRS");

  Tensor<> realXRaij(4, Rvoo, syms, *PirR->wrld, "RealXRaij");
  Tensor<> imagXRaij(4, Rvoo, syms, *PirR->wrld, "ImagXRaij");

  // Allocate and compute PiaR
  int aRStart[] = {No , 0};
  int aREnd[]   = {Np ,NR};
  Tensor<complex> PiaR(PirR->slice(aRStart,aREnd));
  PiaR.set_name("PiaR");

  // Slice the respective parts from PiaR
  int leftPiStart[]  = { 0 ,                                a };
  int leftPiEnd[]    = { Nv, std::min(a+factorsSliceSize, NR) };
  int rightPiStart[] = { 0 ,                                b };
  int rightPiEnd[]   = { Nv, std::min(b+factorsSliceSize, NR) };

  Tensor<complex> leftPiaR (PiaR.slice(leftPiStart  ,  leftPiEnd));
  leftPiaR.set_name("leftPiaR");
  Tensor<complex> rightPiaR(PiaR.slice(rightPiStart , rightPiEnd));
  rightPiaR.set_name("rightPiaR");

  // Split left and right PiaR into real and imaginary parts
  Tensor<> realLeftPiaR(2, leftPiaR.lens, leftPiaR.sym, *leftPiaR.wrld, "RealLeftPiaR");
  Tensor<> imagLeftPiaR(2, leftPiaR.lens, leftPiaR.sym, *leftPiaR.wrld, "ImagRightPiaR");
  fromComplexTensor(leftPiaR, realLeftPiaR, imagLeftPiaR);

  Tensor<> realRightPiaR(2, rightPiaR.lens, rightPiaR.sym, *rightPiaR.wrld, "RealLeftPiaR");
  Tensor<> imagRightPiaR(2, rightPiaR.lens, rightPiaR.sym, *rightPiaR.wrld, "ImagRightPiaR");
  fromComplexTensor(leftPiaR, realLeftPiaR, imagLeftPiaR);

  // Slice the respective parts from LambdaGR
  int leftLambdaStart[]  = { 0  ,                            a };
  int leftLambdaEnd[]    = { NG , std::min(a+factorsSliceSize, NR) };
  Tensor<complex> leftLambdaGR (LambdaGR->slice(leftLambdaStart , leftLambdaEnd));
  leftLambdaGR.set_name("leftLambdaGR");

  int rightLambdaStart[]  = { 0  ,                            b };
  int rightLambdaEnd[]    = { NG , std::min(b+factorsSliceSize, NR) };
  Tensor<complex> rightLambdaGR (LambdaGR->slice(rightLambdaStart , rightLambdaEnd));
  rightLambdaGR.set_name("rightLambdaGR");

  // FIXME: Currently assuming GammaGqr = PirR*PirR*LambdaGR
  //        First Pi not conjugated.
  realXRaij["Rdij"] = (+1.0) * Iabij["cdij"] * realLeftPiaR["cR"];
  imagXRaij["Rdij"] = (-1.0) * Iabij["cdij"] * imagLeftPiaR["cR"];
  Tensor<complex> XRaij(4, Rvoo, syms, *PirR->wrld, "XRaij");
  toComplexTensor(realXRaij, imagXRaij, XRaij);

  Tensor<complex> XRSij(4, RRoo, syms, *PirR->wrld, "XRSij");
  XRSij["RSij"] = XRaij["Rdij"] * rightPiaR["dS"];

  Univar_Function<complex> fConj(&cc4s::conj<complex>);
  Tensor<complex> conjLeftLambdaGR(false, leftLambdaGR);
  conjLeftLambdaGR.set_name("ConjLeftLambdaGR");
  conjLeftLambdaGR.sum(1.0, leftLambdaGR,"GR", 0.0,"GR", fConj);
  VRS["RS"] = conjLeftLambdaGR["GR"] * rightLambdaGR["GS"];

  XRSij["RSij"] = XRSij["RSij"]  * VRS["RS"];

  // Allocate and compute PiiR
  int iRStart[] = {0 , 0};
  int iREnd[]   = {No ,NR};
  Tensor<complex> PiiR(PirR->slice(iRStart,iREnd));
  PiiR.set_name("PiiR");

  // Split PiiR into real and imaginary parts
  Tensor<> realPiiR(2, PiiR.lens, PiiR.sym, *PiiR.wrld, "RealPiiR");
  Tensor<> imagPiiR(2, PiiR.lens, PiiR.sym, *PiiR.wrld, "ImagPiiR");
  fromComplexTensor(PiiR, realPiiR, imagPiiR);

  // Initialize dressedPiaR
  Tensor<complex> dressedPiaR(PiaR);
  dressedPiaR.set_name("dressedPiaR");

  // Split dressedPiaR into real and imaginary parts
  Tensor<> realDressedPiaR(2, dressedPiaR.lens, dressedPiaR.sym, *dressedPiaR.wrld, "RealDressedPiaR");
  Tensor<> imagDressedPiaR(2, dressedPiaR.lens, dressedPiaR.sym, *dressedPiaR.wrld, "ImagDressedPiaR");
  fromComplexTensor(dressedPiaR, realDressedPiaR, imagDressedPiaR);

  // Construct dressedPiaR
  realDressedPiaR["aR"] += (-1.0) * realPiiR["kR"] * (*Tai)["ak"];
  imagDressedPiaR["aR"] += (-1.0) * imagPiiR["kR"] * (*Tai)["ak"];
  toComplexTensor(realDressedPiaR, imagDressedPiaR, dressedPiaR);

  // Slice the respective parts from dressedPiaR
  Tensor<complex> dressedLeftPiaR (dressedPiaR.slice(leftPiStart  ,  leftPiEnd));
  dressedLeftPiaR.set_name("dressedLeftPiaR");
  Tensor<complex> dressedRightPiaR(dressedPiaR.slice(rightPiStart , rightPiEnd));
  dressedRightPiaR.set_name("dressedRightPiaR");

  // Split dressed left PiaR into real and imaginary parts
  Tensor<> dressedRealLeftPiaR(2, dressedLeftPiaR.lens, dressedLeftPiaR.sym, *dressedLeftPiaR.wrld, "dressedRealLeftPiaR");
  Tensor<> dressedImagLeftPiaR(2, dressedLeftPiaR.lens, dressedLeftPiaR.sym, *dressedLeftPiaR.wrld, "dressedImagLeftPiaR");
  fromComplexTensor(dressedLeftPiaR, dressedRealLeftPiaR, dressedImagLeftPiaR);

  XRaij["Rbij"] = XRSij["RSij"]  * dressedRightPiaR["bS"];

  // allocate Tensor for sliced T2 amplitudes
  int vvoo[] = { Nv, Nv, No, No };
  Tensor<> *Fabij(new Tensor<>(4, vvoo, syms, *PirR->wrld, "Fabij"));

  // compute sliced amplitudes
  fromComplexTensor(XRaij, realXRaij, imagXRaij);
  (*Fabij)["abij"]  = realXRaij["Rbij"]  * dressedRealLeftPiaR["aR"];
  (*Fabij)["abij"] += imagXRaij["Rbij"]  * dressedImagLeftPiaR["aR"];

  // return sliced amplitudes
  return Fabij;
}

Tensor<> *ClusterSinglesDoublesAlgorithm::amplitudesFromCoupledCoulombFactors()
{

  // Read the doubles amplitudes Tabij
  Tensor<> *Tabij(&TabijMixer->getNext());
  Tabij->set_name("Tabij");
  Tensor<> *Tai(&TaiMixer->getNext());
  Tai->set_name("Tai");

  // Intermediate tensor Iabij=T2+T1*T1
  Tensor<> Iabij(Tabij);
  Iabij.set_name("Iabij");
  Iabij["abij"] += (*Tai)["ai"] * (*Tai)["bj"];

  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Read the Coulomb Factors PiqR and LambdaGR
  Tensor<complex> *PiqR(getTensorArgument<complex>("FactorOrbitals"));
  PiqR->set_name("PiqR");
  Tensor<complex> *LambdaGR(getTensorArgument<complex>("CoulombFactors"));
  LambdaGR->set_name("LambdaGR");

  int Np=No+Nv;
  int NR(PiqR->lens[1]);
  int Rvoo[] = { NR, Nv, No, No };
  int RRoo[] = { NR, NR, No, No };
  int RR[] = { NR, NR };
  int syms[] = { NS, NS, NS, NS };

  Tensor<complex> VRS(2, RR, syms, *epsi->wrld, "VRS");

  Tensor<> realXRaij(4, Rvoo, syms, *epsi->wrld, "RealXRaij");
  Tensor<> imagXRaij(4, Rvoo, syms, *epsi->wrld, "ImagXRaij");

  // Allocate and compute PiaR
  int aRStart[] = {No , 0};
  int aREnd[]   = {Np ,NR};
  Tensor<complex> PiaR(PiqR->slice(aRStart,aREnd));
  PiaR.set_name("PiaR");

  // Split PiaR into real and imaginary parts
  Tensor<> realPiaR(2, PiaR.lens, PiaR.sym, *PiaR.wrld, "RealPiaR");
  Tensor<> imagPiaR(2, PiaR.lens, PiaR.sym, *PiaR.wrld, "ImagPiaR");
  fromComplexTensor(PiaR, realPiaR, imagPiaR);

  // FIXME: Currently assuming GammaGqr = PiqR*PirR*LambdaGR
  //        First Pi not conjugated.
  realXRaij["Rdij"] = +1.0 * Iabij["cdij"] * realPiaR["cR"];
  imagXRaij["Rdij"] = -1.0 * Iabij["cdij"] * imagPiaR["cR"];
  Tensor<complex> XRaij(4, Rvoo, syms, *epsi->wrld, "XRaij");
  toComplexTensor(realXRaij, imagXRaij, XRaij);

  Tensor<complex> XRSij(4, RRoo, syms, *epsi->wrld, "XRSij");
  XRSij["RSij"] = XRaij["Rdij"] * PiaR["dS"];

  Univar_Function<complex> fConj(&cc4s::conj<complex>);
  Tensor<complex> conjLambdaGR(false, *LambdaGR);
  // conjLambdaGR["GR"] = conj(LambdaGR["GR"])
  conjLambdaGR.set_name("ConjLambdaGR");
  conjLambdaGR.sum(1.0, *LambdaGR,"GR", 0.0,"GR", fConj);
  VRS["RS"] = conjLambdaGR["GR"] * (*LambdaGR)["GS"];

  XRSij["RSij"] = XRSij["RSij"] * VRS["RS"];

  // Allocate and compute PiiR
  int iRStart[] = {0 , 0};
  int iREnd[]   = {No ,NR};
  Tensor<complex> PiiR(PiqR->slice(iRStart,iREnd));
  PiiR.set_name("PiiR");

  // Split PiiR into real and imaginary parts
  Tensor<> realPiiR(2, PiiR.lens, PiiR.sym, *PiiR.wrld, "RealPiiR");
  Tensor<> imagPiiR(2, PiiR.lens, PiiR.sym, *PiiR.wrld, "ImagPiiR");
  fromComplexTensor(PiiR, realPiiR, imagPiiR);

  // Initialize dressedPiaR
  Tensor<complex> dressedPiaR(PiaR);
  dressedPiaR.set_name("dressedPiaR");

  // Split dressedPiaR into real and imaginary parts
  Tensor<> realDressedPiaR(2, dressedPiaR.lens, dressedPiaR.sym, *dressedPiaR.wrld, "RealDressedPiaR");
  Tensor<> imagDressedPiaR(2, dressedPiaR.lens, dressedPiaR.sym, *dressedPiaR.wrld, "ImagDressedPiaR");
  fromComplexTensor(dressedPiaR, realDressedPiaR, imagDressedPiaR);

  // Construct dressedPiaR
  realDressedPiaR["aR"] += (-1.0) * realPiiR["kR"] * (*Tai)["ak"];
  imagDressedPiaR["aR"] += (-1.0) * imagPiiR["kR"] * (*Tai)["ak"];
  toComplexTensor(realDressedPiaR, imagDressedPiaR, dressedPiaR);

  // Contract dressedPiaR with XRSij
  XRaij["Rbij"] = XRSij["RSij"]  * dressedPiaR["bS"];
  fromComplexTensor(XRaij, realXRaij, imagXRaij);

  // allocate Tensor for T2 amplitudes
  int vvoo[] = { Nv, Nv, No, No };
  Tensor<> *Fabij(new Tensor<>(4, vvoo, syms, *PiqR->wrld, "Fabij"));

  // compute sliced amplitudes
  (*Fabij)["abij"] += realXRaij["Rbij"]  * realDressedPiaR["aR"];
  (*Fabij)["abij"] += imagXRaij["Rbij"]  * imagDressedPiaR["aR"];

  // return sliced amplitudes
  return Fabij;
}

DryTensor<> *ClusterSinglesDoublesAlgorithm::drySliceAmplitudesFromCoupledCoulombFactors(int factorsSliceSize)
{
  getTensorArgument<complex,DryTensor<complex>>("FactorOrbitals");
  DryTensor<complex> *LambdaGR(getTensorArgument<complex, 
                               DryTensor<complex>>("CoulombFactors"));
  
  DryTensor<> *epsa(getTensorArgument
                    <double, DryTensor<double>>("ParticleEigenEnergies"));
  DryTensor<> *epsi(getTensorArgument
                    <double, DryTensor<double>>("HoleEigenEnergies"));

  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int NG(LambdaGR->lens[0]);
  int Rvoo[] = { factorsSliceSize,               Nv, No, No };
  int RRoo[] = { factorsSliceSize, factorsSliceSize, No, No };
  int vvoo[] = {               Nv,               Nv, No, No };
  int   vo[] = {               Nv,               No,        };
  int   vR[] = {               Nv, factorsSliceSize         };
  int   RR[] = { factorsSliceSize, factorsSliceSize         };
  int   GR[] = {               NG, factorsSliceSize         };
  int syms[] = {               NS,               NS, NS, NS };

  // Read the doubles amplitudes Tabij
  DryTensor<> Tabij(4, vvoo, syms, SOURCE_LOCATION);
  DryTensor<>   Tai(2,   vo, syms, SOURCE_LOCATION);

  DryTensor<> Iabij(4, vvoo, syms, SOURCE_LOCATION);

  DryTensor<complex> VRS(2, RR, syms, SOURCE_LOCATION);

  DryTensor<> realXRaij(4, Rvoo, syms, SOURCE_LOCATION);
  DryTensor<> imagXRaij(4, Rvoo, syms, SOURCE_LOCATION);

  // Allocate PiaR
  DryTensor<complex> PiaR(2, vR, syms, SOURCE_LOCATION);

  // Slice the respective parts from PiaR
  DryTensor<complex>  leftPiaR(2, vR, syms, SOURCE_LOCATION);
  DryTensor<complex> rightPiaR(2, vR, syms, SOURCE_LOCATION);

  // Split left and right PiaR into real and imaginary parts
  DryTensor<>  realLeftPiaR(2, vR, syms, SOURCE_LOCATION);
  DryTensor<>  imagLeftPiaR(2, vR, syms, SOURCE_LOCATION);
  DryTensor<> realRightPiaR(2, vR, syms, SOURCE_LOCATION);
  DryTensor<> imagRightPiaR(2, vR, syms, SOURCE_LOCATION);

  // Slice the respective parts from LambdaGR
  DryTensor<complex>  leftLambdaGR(2, GR, syms, SOURCE_LOCATION);
  DryTensor<complex> rightLambdaGR(2, GR, syms, SOURCE_LOCATION);

  // Allocate dressed PiaR
  DryTensor<complex> dressedPiaR(2, vR, syms, SOURCE_LOCATION);

  // Slice the respective parts from dressed PiaR
  DryTensor<complex>  leftDressedPiaR(2, vR, syms, SOURCE_LOCATION);
  DryTensor<complex> rightDressedPiaR(2, vR, syms, SOURCE_LOCATION);

  // Split left and right PiaR into real and imaginary parts
  DryTensor<>  realLeftDressedPiaR(2, vR, syms, SOURCE_LOCATION);
  DryTensor<>  imagLeftDressedPiaR(2, vR, syms, SOURCE_LOCATION);

  // FIXME: Currently assuming GammaGqr = PirR*PirR*LambdaGR
  //        First Pi not conjugated.
  DryTensor<complex> XRaij(4, Rvoo, syms, SOURCE_LOCATION);

  DryTensor<complex> XRSij(4, RRoo, syms, SOURCE_LOCATION);

  DryTensor<complex> conjLeftLambdaGR(leftLambdaGR, SOURCE_LOCATION);

  // allocate Tensor for sliced T2 amplitudes
  DryTensor<> *Fabij(new DryTensor<>(4, vvoo, syms, SOURCE_LOCATION));

  // return sliced amplitudes
  return Fabij;
}
