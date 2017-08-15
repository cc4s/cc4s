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
  int excitationLevel(D.order/2);
  for (int p(0); p < excitationLevel; ++p) {
    char epsiIndex[2] = {indices[excitationLevel+p], 0};
    D[indices.c_str()] += Fepsi[epsiIndex];
    char epsaIndex[2] = {indices[p], 0};
    D[indices.c_str()] -= Fepsa[epsaIndex];
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


Tensor<double> *ClusterSinglesDoublesAlgorithm::sliceCoupledCoulombIntegrals(
  Mixer<double> *TaiMixer,
  int a, int b, int integralsSliceSize
) {
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

Tensor<complex> *ClusterSinglesDoublesAlgorithm::sliceCoupledCoulombIntegrals(
  Mixer<complex> *TaiMixer,
  int a, int b, int integralsSliceSize
) {
  // Read the amplitudes Tai
  Tensor<complex> *Tai(&TaiMixer->getNext());
  Tai->set_name("Tai");

  // Read the Coulomb vertex GammaGqr
  Tensor<complex> *GammaGqr( getTensorArgument<complex>("CoulombVertex"));
  GammaGqr->set_name("GammaGqr");

  // Compute No,Nv,NG,Np
  int No(Tai->lens[1]);
  int Nv(Tai->lens[0]);
  int NG(GammaGqr->lens[0]);
  int Np(No+Nv);

  int aStart(Np-Nv), aEnd(Np);
  int iStart(0), iEnd(No);
  int GiaStart[] = {0, iStart,aStart};
  int GiaEnd[]   = {NG,iEnd,  aEnd};
  int GaiStart[] = {0, aStart,iStart};
  int GaiEnd[]   = {NG,aEnd,  iEnd};
  int GabStart[] = {0, aStart,aStart};
  int GabEnd[]   = {NG,aEnd,  aEnd};
  auto GammaGia( new Tensor<complex>(GammaGqr->slice(GiaStart, GiaEnd)) );
  auto GammaGai( new Tensor<complex>(GammaGqr->slice(GaiStart, GaiEnd)) );
  auto GammaGab( new Tensor<complex>(GammaGqr->slice(GabStart, GabEnd)) );

  Univar_Function<complex> fConj(conj<complex>);

  Tensor<complex> conjTransposeGammaGia(false, *GammaGia);
  conjTransposeGammaGia.sum(1.0,*GammaGai,"Gai", 0.0,"Gia", fConj);
  Tensor<complex> conjTransposeGammaGab(false, *GammaGab);
  conjTransposeGammaGab.sum(1.0,*GammaGab,"Gba", 0.0,"Gab", fConj);

  // Construct dressed Coulomb vertex GammaGab
  Tensor<complex> DressedGammaGab(GammaGab);
  DressedGammaGab.set_name("DressedGammaGab");
  DressedGammaGab["Gab"] += (-1.0) * (*GammaGia)["Gkb"] * (*Tai)["ak"];

  Tensor<complex> conjTransposeDressedGammaGab(conjTransposeGammaGab);
  conjTransposeDressedGammaGab.set_name("conjTransposeDressedGammaGab");
  conjTransposeDressedGammaGab["Gab"] += (-1.0) * conjTransposeGammaGia["Gkb"] * (*Tai)["ak"];

  // Slice the respective parts from the dressed Coulomb vertex GammaGab
  int leftGammaStart[] = { 0, a, 0 };
  int leftGammaEnd[] = { NG, std::min(a+integralsSliceSize, Nv), Nv };
  int rightGammaStart[] = { 0, b, 0 };
  int rightGammaEnd[] = { NG, std::min(b+integralsSliceSize, Nv), Nv };

  Tensor<complex> leftGamma(conjTransposeDressedGammaGab.slice(leftGammaStart, leftGammaEnd));
  Tensor<complex> rightGamma(DressedGammaGab.slice(rightGammaStart, rightGammaEnd));

  // Allocate sliced Coulomb integrals
  int lens[] = {
    leftGamma.lens[1], rightGamma.lens[1], leftGamma.lens[2], rightGamma.lens[2]
  };
  int syms[] = { NS, NS, NS, NS };
  Tensor<complex> *Vxycd(new Tensor<complex>(4, lens, syms, *GammaGqr->wrld, "Vxycd"));

  // Contract left and right slices of the dressed Coulomb vertices
  (*Vxycd)["xycd"]  = leftGamma["Gxc"] * rightGamma["Gyd"];
  return Vxycd;
}


Tensor<double> *
  ClusterSinglesDoublesAlgorithm::sliceAmplitudesFromCoupledCoulombFactors
(
  Mixer<double> *TaiMixer, Mixer<double> *TabijMixer,
  int a, int b, int factorsSliceSize
) {
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

  // TODO: specify how the vertex should be computed
  // assuming GammaGqr = PirR*PirR*LambdaGR (first Pi not conjugated)
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

Tensor<complex> *
  ClusterSinglesDoublesAlgorithm::sliceAmplitudesFromCoupledCoulombFactors
(
  Mixer<complex> *TaiMixer, Mixer<complex> *TabijMixer,
  int a, int b, int factorsSliceSize
) {
  Tensor<complex> *PirR(getTensorArgument<complex>("FactorOrbitals"));
  PirR->set_name("PirR");
  Tensor<complex> *LambdaGR(getTensorArgument<complex>("CoulombFactors"));
  LambdaGR->set_name("LambdaGR");

  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));

  // Read the doubles amplitudes Tabij
  Tensor<complex> *Tabij(&TabijMixer->getNext());
  Tabij->set_name("Tabij");
  Tensor<complex> *Tai(&TaiMixer->getNext());
  Tai->set_name("Tai");

  // Intermediate tensor Iabij=T2+T1*T1
  Tensor<complex> Iabij(Tabij);
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

  Univar_Function<complex> fConj(&cc4s::conj<complex>);

  Tensor<complex> VRS(2, RR, syms, *PirR->wrld, "VRS");

  Tensor<complex> XRaij(4, Rvoo, syms, *PirR->wrld, "XRaij");

  // Allocate and compute PiaR
  int aRStart[] = {No , 0};
  int aREnd[]   = {Np ,NR};
  Tensor<complex> PiaR(PirR->slice(aRStart,aREnd));
  PiaR.set_name("PiaR");
  Tensor<complex> conjPiaR(false, PiaR);
  conjPiaR.set_name("ConjPiaR");
  conjPiaR.sum(1.0, PiaR,"aR", 0.0,"aR", fConj);

  // Slice the respective parts from PiaR
  int leftPiStart[]  = { 0 ,                                a };
  int leftPiEnd[]    = { Nv, std::min(a+factorsSliceSize, NR) };
  int rightPiStart[] = { 0 ,                                b };
  int rightPiEnd[]   = { Nv, std::min(b+factorsSliceSize, NR) };

  Tensor<complex> leftPiaR (PiaR.slice(leftPiStart  ,  leftPiEnd));
  leftPiaR.set_name("leftPiaR");
  Tensor<complex> rightPiaR(PiaR.slice(rightPiStart , rightPiEnd));
  rightPiaR.set_name("rightPiaR");
  
  // Slice the respective parts from LambdaGR
  int leftLambdaStart[]  = { 0  ,                                a };
  int leftLambdaEnd[]    = { NG , std::min(a+factorsSliceSize, NR) };
  Tensor<complex> leftLambdaGR (LambdaGR->slice(leftLambdaStart , leftLambdaEnd));
  leftLambdaGR.set_name("leftLambdaGR");

  Tensor<complex> conjLeftLambdaGR(false, leftLambdaGR);
  conjLeftLambdaGR.set_name("ConjLeftLambdaGR");
  conjLeftLambdaGR.sum(1.0, leftLambdaGR,"GR", 0.0,"GR", fConj);

  int rightLambdaStart[]  = { 0  ,                                b };
  int rightLambdaEnd[]    = { NG , std::min(b+factorsSliceSize, NR) };
  Tensor<complex> rightLambdaGR (LambdaGR->slice(rightLambdaStart , rightLambdaEnd));
  rightLambdaGR.set_name("rightLambdaGR");

  // TODO: specify how the vertex should be computed
  // assuming GammaGqr = (PiqR*)*(PirR)*(LambdaGR) (first Pi conjugated)
  XRaij["Rdij"] = (+1.0) * Iabij["cdij"] * leftPiaR["cR"];

  Tensor<complex> XRSij(4, RRoo, syms, *PirR->wrld, "XRSij");
  XRSij["RSij"] = XRaij["Rdij"] * rightPiaR["dS"];

  VRS["RS"] = conjLeftLambdaGR["GR"] * rightLambdaGR["GS"];

  XRSij["RSij"] = XRSij["RSij"]  * VRS["RS"];

  // Allocate and compute PiiR
  int iRStart[] = {0 , 0};
  int iREnd[]   = {No ,NR};
  Tensor<complex> PiiR(PirR->slice(iRStart,iREnd));
  PiiR.set_name("PiiR");
  Tensor<complex> conjPiiR(false, PiiR);
  conjPiiR.set_name("ConjPiiR");
  conjPiiR.sum(1.0, PiiR,"iR", 0.0,"iR", fConj);

  // Initialize dressedPiaR
  Tensor<complex> dressedPiaR(conjPiaR);
  dressedPiaR.set_name("dressedPiaR");

  // Construct dressedPiaR
  dressedPiaR["aR"] += (-1.0) * conjPiiR["kR"] * (*Tai)["ak"];

  // Slice the respective parts from dressedPiaR
  Tensor<complex> dressedLeftPiaR (dressedPiaR.slice(leftPiStart  ,  leftPiEnd));
  dressedLeftPiaR.set_name("dressedLeftPiaR");
  Tensor<complex> dressedRightPiaR(dressedPiaR.slice(rightPiStart , rightPiEnd));
  dressedRightPiaR.set_name("dressedRightPiaR");

  XRaij["Rbij"] = XRSij["RSij"]  * dressedRightPiaR["bS"];

  // allocate Tensor for sliced T2 amplitudes
  int vvoo[] = { Nv, Nv, No, No };
  Tensor<complex> *Fabij(new Tensor<complex>(4, vvoo, syms, *PirR->wrld, "Fabij"));

  // compute sliced amplitudes
  (*Fabij)["abij"]  = XRaij["Rbij"]  * dressedLeftPiaR["aR"];

  // return sliced amplitudes
  return Fabij;
}


Tensor<double> *
  ClusterSinglesDoublesAlgorithm::amplitudesFromCoupledCoulombFactors
(
  Mixer<double> *TaiMixer, Mixer<double> *TabijMixer
) {
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

Tensor<complex> *
  ClusterSinglesDoublesAlgorithm::amplitudesFromCoupledCoulombFactors
(
  Mixer<complex> *TaiMixer, Mixer<complex> *TabijMixer
) {
  // FIXME: implement me for complex amplitudes...
  return nullptr;
}


template <typename F>
void ClusterSinglesDoublesAlgorithm::sliceIntoResiduum(
  Tensor<F> &Rxyij, int a, int b, Tensor<F> &Rabij
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
    Tensor<F> Ryxji(4, srcEnd, Rxyij.sym, *Rxyij.wrld, "Ryxji");
    Ryxji["yxji"] = Rxyij["xyij"];
    // Add Ryxij to Rabij
    Rabij.slice(dstStart,dstEnd,1.0, Ryxji,srcStart,srcEnd,1.0);
  }
}

// instantiate:
template
void ClusterSinglesDoublesAlgorithm::sliceIntoResiduum(
  Tensor<double> &Rxyij, int a, int b, Tensor<double> &Rabij
);
template
void ClusterSinglesDoublesAlgorithm::sliceIntoResiduum(
  Tensor<complex> &Rxyij, int a, int b, Tensor<complex> &Rabij
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

