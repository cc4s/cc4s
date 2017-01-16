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
    throw new EXCEPTION(stringStream.str());
  }

  {
    // Allocate the doubles amplitudes and append it to the mixer
    int No(Vabij->lens[2]);
    int Nv(Vabij->lens[0]);
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
  if (isArgumentGiven(amplitudesName.str())) {
    allocatedTensorArgument(
      amplitudesName.str(), new Tensor<>(TabijMixer->getNext())
    );
  }

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
    throw new EXCEPTION(stringStream.str());
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
    DryTensor<> Tabij(4, vvoo, syms, SOURCE_LOCATION);
    allocatedTensorArgument(
      amplitudesName.str(), new DryTensor<>(Tabij, SOURCE_LOCATION)
    );
  }

  // Allocate the energy e
  DryScalar<> energy(SOURCE_LOCATION);

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
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));

  // Build Dabij
  Tensor<> Dabij(false, Rabij);
  Dabij.set_name("Dabij");
  Dabij["abij"]  = (*epsi)["i"];
  Dabij["abij"] += (*epsi)["j"];
  Dabij["abij"] -= (*epsa)["a"];
  Dabij["abij"] -= (*epsa)["b"];

  // TODO:
  // levelshifting can be implemented here

  // Divide Rabij/Dabij to get Tabij
  Bivar_Function<> fDivide(&divide<double>);
  Rabij.contract(1.0, Rabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);

/*
  // TODO: test why bivariate transform doesn't work
  Matrix<> Dai(epsa->lens[0], epsi->lens[0], *epsi->wrld);
  Dai.set_name("Dai");
  Dai["ai"]  = (*epsi)["i"];
  Dai["ai"] -= (*epsa)["a"];
  Transform<double, double, double>(
    std::function<void(double, double, double &)>(
      [](double Dai, double Dbj, double &R) {
        R /= Dai + Dbj;
      }
    )
  ) (
    Dai["ai"], Dai["bj"], Rabij["abij"]
  );
*/
}

void ClusterDoublesAlgorithm::dryDoublesAmplitudesFromResiduum(
  cc4s::DryTensor<> &Rabij
) {
  // Build Dabij
  DryTensor<> Dabij(Rabij, SOURCE_LOCATION);
}


Tensor<> *ClusterDoublesAlgorithm::sliceCoulombIntegrals(int a, int b, int integralsSliceSize) {
  Tensor<complex> *GammaGqr(getTensorArgument<complex>("CoulombVertex"));
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(No+Nv);
  int NG(GammaGqr->lens[0]);
  
  // Slice the respective parts from the Coulomb vertex
  int leftGammaStart[] = { 0, No+a, No };
  int leftGammaEnd[] = { NG, std::min(No+a+integralsSliceSize, Np), Np };
  int rightGammaStart[] = { 0, No+b, No };
  int rightGammaEnd[] = { NG, std::min(No+b+integralsSliceSize, Np), Np };
  // FIXME: replace copy constructor calls by allocation then T.slice(...) call
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

DryTensor<> *ClusterDoublesAlgorithm::drySliceCoulombIntegrals(int integralsSliceSize) {
  // Read the Coulomb vertex GammaGqr
  DryTensor<complex> *GammaGqr(getTensorArgument<complex, 
                               DryTensor<complex>>("CoulombVertex"));
  
  // Read the Particle/Hole Eigenenergies
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies")
  );

  int Nv(epsa->lens[0]);
  int NG(GammaGqr->lens[0]);
  
  // Slice the respective parts from the Coulomb vertex
  int syms[] = { NS, NS, NS, NS };
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

void ClusterDoublesAlgorithm::printEnergyFromResiduum(CTF::Tensor<> &Rabij)
{
  // Read the Coulomb Integrals Vabij required for the energy
  Tensor<> *Vabij(getTensorArgument<>("PPHHCoulombIntegrals"));

  Tensor<> Eabij(Rabij);
  Eabij.set_name("Eabij");

  doublesAmplitudesFromResiduum(Eabij);

  // Allocate the energy e
  Scalar<> energy(*Vabij->wrld);
  energy.set_name("energy");
  double e(0), dire, exce;

  std::string abbreviation(getAbbreviation());
  std::transform(abbreviation.begin(), abbreviation.end(), 
                 abbreviation.begin(), ::toupper);

  // Direct term
  energy[""] = 2.0 * Eabij["abij"] * (*Vabij)["abij"];
  dire  = energy.get_val();
  // Exchange term
  energy[""] = Eabij["abji"] * (*Vabij)["abij"];
  exce  = -1.0 * energy.get_val();
  // Total energy
  e = dire + exce;
  LOG(2, abbreviation) << "Energy=" << e << std::endl;
}

Tensor<> *ClusterDoublesAlgorithm::sliceAmplitudesFromCoulombFactors(int a,
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
  int leftPiStart[]  = { 0 ,                            a };
  int leftPiEnd[]    = { Nv, std::min(a+factorsSliceSize, NR) };
  int rightPiStart[] = { 0 ,                            b };
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
  realXRaij["Rdij"] = (+1.0) * (*Tabij)["cdij"] * realLeftPiaR["cR"];
  imagXRaij["Rdij"] = (-1.0) * (*Tabij)["cdij"] * imagLeftPiaR["cR"];
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
  XRaij["Rbij"] = XRSij["RSij"]  * rightPiaR["bS"];

  // allocate Tensor for sliced T2 amplitudes
  int vvoo[] = { Nv, Nv, No, No };
  Tensor<> *Xabij(new Tensor<>(4, vvoo, syms, *PirR->wrld, "Xabij"));

  // compute sliced amplitudes
  fromComplexTensor(XRaij, realXRaij, imagXRaij);
  (*Xabij)["abij"] += realXRaij["Rbij"]  * realLeftPiaR["aR"];
  (*Xabij)["abij"] += imagXRaij["Rbij"]  * imagLeftPiaR["aR"];

  // return sliced amplitudes
  return Xabij;
}

DryTensor<> *ClusterDoublesAlgorithm::drySliceAmplitudesFromCoulombFactors(int factorsSliceSize)
{
  getTensorArgument<complex,DryTensor<complex>>("FactorOrbitals");
  DryTensor<complex> *LambdaGR(getTensorArgument<complex, 
                               DryTensor<complex>>("CoulombFactors"));
  
  // Read the Particle/Hole Eigenenergies
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies")
  );
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies")
  );

  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int NG(LambdaGR->lens[0]);
  int Rvoo[] = { factorsSliceSize,               Nv, No, No };
  int RRoo[] = { factorsSliceSize, factorsSliceSize, No, No };
  int vvoo[] = {               Nv,               Nv, No, No };
  int   vR[] = {               Nv, factorsSliceSize         };
  int   RR[] = { factorsSliceSize, factorsSliceSize         };
  int   GR[] = {               NG, factorsSliceSize         };
  int syms[] = {               NS,               NS, NS, NS };

  // Read the doubles amplitudes Tabij
  DryTensor<> Tabij(4, vvoo, syms, SOURCE_LOCATION);

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


Tensor<> *ClusterDoublesAlgorithm::sliceAmplitudesFromCoulombFactorsTcc(int a,
                                                                        int factorsSliceSize)
{
  Tensor<complex> *PirR(getTensorArgument<complex>("FactorOrbitals"));
  PirR->set_name("PirR");
  Tensor<complex> *LambdaGR(getTensorArgument<complex>("CoulombFactors"));
  LambdaGR->set_name("LambdaGR");

  Tensor<complex> *GammaGqr(getTensorArgument<complex>("CoulombVertex"));

  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));

  // Read the doubles amplitudes Tabij
  Tensor<> *Tabij(&TabijMixer->getNext());
  Tabij->set_name("Tabij");

  // Define complex conjugation function
  Univar_Function<complex> fConj(&cc4s::conj<complex>);

  // Define tcc environment
  auto machineTensorFactory(CtfMachineTensor<complex>::Factory::create());
  auto tcc(tcc::Tcc<complex>::create(machineTensorFactory));

  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(PirR->lens[0]);
  int NR(PirR->lens[1]);
  int NG(LambdaGR->lens[0]);
  int Rx(std::min(factorsSliceSize, NR-a));
  int syms[] = { NS, NS, NS, NS };

  // Allocate and compute GammaGab from GammaGqr
  int iStart(0), iEnd(No);
  int aStart(Np-Nv), aEnd(Np);
  int GabStart[] = {0 ,aStart,aStart};
  int GabEnd[]   = {NG,aEnd,  aEnd  };
  Tensor<complex> GammaGab(GammaGqr->slice(GabStart,GabEnd));

  // Create tcc tensors from ctf tensors (conj)GammaGab
  auto tccGammaGab( tcc->createTensor(CtfMachineTensor<complex>::create(GammaGab)));

  // Allocate and compute PiaR
  int aRStart[] = {No , 0};
  int aREnd[]   = {Np ,NR};
  Tensor<complex> PiaR(PirR->slice(aRStart,aREnd));
  PiaR.set_name("PiaR");

  // Slice the respective parts from PiaR
  int leftPiStart[]  = { 0 ,                            a };
  int leftPiEnd[]    = { Nv, std::min(a+factorsSliceSize, NR) };

  Tensor<complex> leftPiaR (PiaR.slice(leftPiStart  ,  leftPiEnd));
  leftPiaR.set_name("leftPiaR");

  Tensor<complex> conjLeftPiaR(false, leftPiaR);
  conjLeftPiaR.set_name("ConjLeftPiaR");
  conjLeftPiaR.sum(1.0, leftPiaR,"aR", 0.0,"aR", fConj);

  // Create tcc tensors from ctf tensors right(left)PiaR
  auto tccConjLeftPiaR( tcc->createTensor(CtfMachineTensor<complex>::create(conjLeftPiaR)) );

  // Slice the respective parts from LambdaGR
  int leftLambdaStart[]  = { 0  ,                            a };
  int leftLambdaEnd[]    = { NG , std::min(a+factorsSliceSize, NR) };
  Tensor<complex> leftLambdaGR (LambdaGR->slice(leftLambdaStart , leftLambdaEnd));
  leftLambdaGR.set_name("leftLambdaGR");

  Tensor<complex> conjLeftLambdaGR(false, leftLambdaGR);
  conjLeftLambdaGR.set_name("ConjLeftLambdaGR");
  conjLeftLambdaGR.sum(1.0, leftLambdaGR,"GR", 0.0,"GR", fConj);

  // Create tcc tensors from ctf tensors right(left)LambdaGR
  auto tccConjLeftLambdaGR( tcc->createTensor(CtfMachineTensor<complex>::create(conjLeftLambdaGR)));

  // Create complex amplitudes tensor
  Tensor<complex> cTabij(4, Tabij->lens, syms, *Tabij->wrld, "cTabij");
  toComplexTensor(*Tabij, cTabij);

  // Create tcc tensor from ctf tensor Tabij
  auto tccCTabij( tcc->createTensor(CtfMachineTensor<complex>::create(cTabij)));

  // Create complex amplitudes tensor
  Tensor<complex> cXabij(false, cTabij);

  // Create tcc tensor from ctf tensor Tabij
  auto tccCXabij( tcc->createTensor(CtfMachineTensor<complex>::create(cXabij)));

  // Do the tcc contraction
  // compile
  auto operation(
    tcc->compile(
      (*tccCXabij)["abij"] <<=
        (*tccCTabij)["cdij"] *
        (*tccGammaGab)["Gbd"] *
        (*tccConjLeftPiaR)["cR"] * (*tccConjLeftPiaR)["aR"] *
        (*tccConjLeftLambdaGR)["GR"]
    )
  );
  // and execute
  operation->execute();

  // allocate Tensor for sliced T2 amplitudes
  int vvoo[] = { Nv, Nv, No, No };
  Tensor<> *Xabij(new Tensor<>(4, vvoo, syms, *PirR->wrld, "Xabij"));
  // for now: duplicate result
  auto implementationXabij(
    std::dynamic_pointer_cast<CtfMachineTensor<complex>>(tccCXabij->getMachineTensor())
  );
  fromComplexTensor(implementationXabij->tensor, *Xabij);

  // return sliced amplitudes
  return Xabij;
}
