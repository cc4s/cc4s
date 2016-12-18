#include <algorithms/LaplaceMp2Energy.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/Tcc.hpp>
#include <tcc/DryTensor.hpp>
#include <util/CtfMachineTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <Options.hpp>

using namespace cc4s;
using namespace CTF;
using std::shared_ptr;
using std::make_shared;

ALGORITHM_REGISTRAR_DEFINITION(LaplaceMp2Energy);

LaplaceMp2Energy::LaplaceMp2Energy(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

LaplaceMp2Energy::~LaplaceMp2Energy() {
  delete GhRSn; delete GpRSn;
  delete VRS; delete wn;
}

void LaplaceMp2Energy::run() {
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<> *Taun(getTensorArgument("ImaginaryTimePoints"));

  // TODO: slice over imaginary time samples to conserve memory
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Nn(Taun->lens[0]);
  int on[] = { No, Nn };
  int vn[] = { Nv, Nn };
  int syms[] = { NS, NS, NS, NS };
  Tensor<> Pan(2, vn, syms, *epsi->wrld, "Pan");
  Tensor<> Hin(2, on, syms, *epsi->wrld, "Hin");

  double mu(getRealArgument("ChemicalPotential"));
  Pan["an"]  = (*epsa)["a"];
  Pan["an"] -= mu;
  Pan["an"] *= (*Taun)["n"];
  CTF::Transform<double>(
    std::function<void(double &)>([](double &c){ c = std::exp(-c); })
  ) (
    Pan["an"]
  );

  Hin["in"]  = (*epsi)["i"];
  Hin["in"] -= mu;
  Hin["in"] *= (*Taun)["n"];
  CTF::Transform<double>(
    std::function<void(double &)>([](double &c){ c = std::exp(+c); })
  ) (
    Hin["in"]
  );

  Tensor<complex> CPan(2, vn, syms, *epsi->wrld, "CPan");
  toComplexTensor(Pan, CPan);

  Tensor<complex> CHin(2, on, syms, *epsi->wrld, "CHin");
  toComplexTensor(Hin, CHin);

  Tensor<complex> *PirR(getTensorArgument<complex>("FactorOrbitals"));
  PirR->set_name("PirR");
  Tensor<complex> *LambdaGR(getTensorArgument<complex>("CoulombFactors"));
  LambdaGR->set_name("LambdaGR"); 

  int Np(PirR->lens[0]);
  int NR(PirR->lens[1]);
  int RR[]  =     { NR, NR };
  int RRn[] = { NR, NR, Nn };

  {
    // check whether factor orbitals are real
    Tensor<> realPirR(2, PirR->lens, PirR->sym, *PirR->wrld, "realPirR");
    Tensor<> imagPirR(2, PirR->lens, PirR->sym, *PirR->wrld, "imagPirR");
    fromComplexTensor(*PirR, realPirR, imagPirR);
    double imagNorm(imagPirR.norm2());
    if (imagNorm > 1e-8) {
      throw new EXCEPTION(
        "real factor orbitals expecteed. Use (realFactorOrbitals 1) in decomposition."
      );
    }
  }

  VRS =   new Tensor<complex>(2, RR,  syms, *PirR->wrld,   "VRS");
  GpRSn = new Tensor<complex>(3, RRn, syms, *PirR->wrld, "GpRSn");
  GhRSn = new Tensor<complex>(3, RRn, syms, *PirR->wrld, "GhRSn");

  // Allocate and compute PiaR
  int aRStart[] = {No, 0};
  int aREnd[]   = {Np,NR};
  Tensor<complex> PiaR(PirR->slice(aRStart,aREnd));
  PiaR.set_name("PiaR");

  // Allocate and compute PiiR
  int iRStart[] = {0,  0};
  int iREnd[]   = {No,NR};
  Tensor<complex> PiiR(PirR->slice(iRStart,iREnd));
  PiiR.set_name("PiiR");

  Univar_Function<complex> fConj(&cc4s::conj<complex>);
  Tensor<complex> conjLambdaGR(false, LambdaGR);
  conjLambdaGR.set_name("ConjLambdaGR");
  conjLambdaGR.sum(1.0, *LambdaGR,"GR", 0.0,"GR", fConj);
  (*VRS)["RS"] = conjLambdaGR["GR"] * (*LambdaGR)["GS"];
  // FIXME: imaginary part manually discarded, should be zero
//  VRS->sum(-0.5, *VRS,"RS", 0.5,"RS", fConj);
  LOG(1, "MP2") << "Coulomb propagator set up" << std::endl;

  // NOTE: the hole propagator sign is considered for the entire term, not here
  (*GhRSn)["RSn"] = PiiR["iS"] * PiiR["iR"] * CHin["in"];
  LOG(1, "MP2") << "Hole propagator set up" << std::endl;

  (*GpRSn)["RSn"] = PiaR["aS"] * PiaR["aR"] * CPan["an"];
  LOG(1, "MP2") << "Particle propagator set up" << std::endl;


  // get numerical weights
  Tensor<> *realwn(getTensorArgument("ImaginaryTimeWeights"));
  wn = new Tensor<complex>(1, &Nn, syms, *PirR->wrld, "wn");
  toComplexTensor(*realwn, *wn);


  calculateAnalytically();
  double direct(calculateDirectTerm());
  double exchange(calculateExchangeTerm());
  double energy(direct + exchange);

//  LOG(0, "MP2") << "e=" << e << std::endl;
//  LOG(1, "MP2") << "MP2d=" << dire << std::endl;
//  LOG(1, "MP2") << "MP2dLaplace=" << calculateDirectTerm() << std::endl;
//  LOG(1, "MP2") << "MP2dAnalytical=" << calculateDirectTermAnalytically() << std::endl;
//  LOG(1, "MP2") << "MP2x=" << exce << std::endl;

  setRealArgument("Mp2Energy", energy);
}

void LaplaceMp2Energy::dryRun() {
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies")
  );

  DryTensor<> *Tn(
    getTensorArgument<double, DryTensor<double>>("ImaginaryTimePoints")
  );
  DryTensor<> *Wn(
    getTensorArgument<double, DryTensor<double>>("ImaginaryTimeWeights")
  );

  DryTensor<complex> *PirR(
    getTensorArgument<complex,DryTensor<complex>>("FactorOrbitals")
  );
  DryTensor<complex> *LambdaGR(
    getTensorArgument<complex, DryTensor<complex>>("CoulombFactors")
  );

  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Nn(Tn->lens[0]);
  int on[] = { No, Nn };
  int vn[] = { Nv, Nn };
  int syms[] = { NS, NS, NS, NS };
  DryTensor<> Pan(2, vn, syms, SOURCE_LOCATION);
  DryTensor<> Hin(2, on, syms, SOURCE_LOCATION);

//  int Np(PirR->lens[0]);
  int NR(PirR->lens[1]);
  int RR[] =   { NR, NR };
  int RRn[] =  { NR, NR, Nn };

  DryTensor<complex> VRS  (2, RR,  syms, SOURCE_LOCATION);
  DryTensor<complex> GpRSn(3, RRn, syms, SOURCE_LOCATION);
  DryTensor<complex> GhRSn(3, RRn, syms, SOURCE_LOCATION);

  // Allocate and compute PiaR
  int vR[] =   { Nv, NR };
  DryTensor<complex> PiaR(2, vR,  syms, SOURCE_LOCATION);

  // Allocate and compute PiiR
  int oR[] =   { No, NR };
  DryTensor<complex> PiiR(2, oR,  syms, SOURCE_LOCATION);

  DryTensor<complex> conjLambdaGR(*LambdaGR);

  DryTensor<complex> conjPiiR(PiiR);

  DryTensor<complex> conjPiaR(PiaR);
  Wn->use();

  DryScalar<> energy();
}

double LaplaceMp2Energy::calculateDirectTerm() {
  Scalar<complex> energy(*Cc4s::world);
  Tensor<complex> ChiVRSn(false, *GpRSn);
  ChiVRSn["RSn"] = (*GpRSn)["RSn"] * (*GhRSn)["SRn"];
  ChiVRSn["RTn"] = ChiVRSn["RSn"] * (*VRS)["ST"];
  Tensor<complex> conjChiVRSn(false, ChiVRSn);
  // second part of the diagram is the conjugate of the first
  Univar_Function<complex> fConj(&conj<complex>);
  conjChiVRSn.sum(1.0, ChiVRSn,"RTn", 0.0,"RTn", fConj);
  // Imaginary time propagators give 1/(eps_a+eps_b-eps_i-eps_j), so negate it
  energy[""] = -2.0 * (*wn)["n"] * ChiVRSn["RTn"] * conjChiVRSn["TRn"];
  complex direct(energy.get_val());
  LOG(1, "MP2") << "Direct energy computed from Laplace = " << direct << std::endl;
  return -2.0 * std::real(energy.get_val());
}

double LaplaceMp2Energy::calculateExchangeTerm() {
  auto machineTensorFactory(CtfMachineTensor<complex>::Factory::create());
  auto tcc(tcc::Tcc<complex>::create(machineTensorFactory));
  auto Gp(tcc->createTensor(CtfMachineTensor<complex>::create(*GpRSn)));
  auto Gh(tcc->createTensor(CtfMachineTensor<complex>::create(*GhRSn)));
  auto V(tcc->createTensor(CtfMachineTensor<complex>::create(*VRS)));
  auto w(tcc->createTensor(CtfMachineTensor<complex>::create(*wn)));
  auto energy(tcc->createTensor(std::vector<int>(), "energy"));
  tcc->compile(
    (*energy)[""] <<=
      (*w)["n"] *
      (*V)["RS"] * (*V)["TU"] *
      (*Gp)["RTn"] * (*Gh)["TSn"] * (*Gp)["SUn"] * (*Gh)["URn"]
  )->execute();
//  energy[""] =  VRS["RS"] * GpRSn["SUn"] * GhRSn["TSn"] * VRS["TU"] *
//    GpRSn["RTn"] * GhRSn["URn"] * (*Wn)["n"];
//  LOG(1, "MP2") << "Exchange energy computed" << std::endl;
//  exce = -1.0 * energy.get_val();

  CTF::Scalar<complex> ctfEnergy;
  ctfEnergy[""] = std::dynamic_pointer_cast<CtfMachineTensor<complex>>(
    energy->getMachineTensor()
  )->tensor[""];
  return std::real(ctfEnergy.get_val());
}

double LaplaceMp2Energy::calculateAnalytically() {
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<complex> *GammaFqr(getTensorArgument<complex>("CoulombVertex"));
  int NF(GammaFqr->lens[0]);
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(No+Nv);
  int FaiStart[] = {0, No, 0};
  int FaiEnd[]   = {NF,Np,No};
  Tensor<complex> GammaFai(GammaFqr->slice(FaiStart,FaiEnd));
  GammaFai.set_name("GammaFai");
  Tensor<complex> conjGammaFai(false, GammaFai);
  Univar_Function<complex> fConj(&conj<complex>);
  conjGammaFai.sum(1.0, GammaFai,"Fai", 0.0,"Fai", fConj);

  int lens[] = { Nv,Nv,No,No };
  int syms[] = { NS,NS,NS,NS };
  Tensor<complex> cVabij(4, lens, syms, *GammaFqr->wrld, "cVabij");
  cVabij["abij"] = conjGammaFai["Fai"] * GammaFai["Fbj"];

  Tensor<complex> cDabij(false, cVabij);
  Tensor<complex> cepsi(1, &No, syms, *epsi->wrld, "cepsi");
  toComplexTensor(*epsi, cepsi);
  Tensor<complex> cepsa(1, &Nv, syms, *epsa->wrld, "cepsi");
  toComplexTensor(*epsa, cepsa);
  cDabij["abij"] =  cepsi["i"];
  cDabij["abij"] += cepsi["j"];
  cDabij["abij"] -= cepsa["a"];
  cDabij["abij"] -= cepsa["b"];
  Bivar_Function<complex> fDivide(&divide<complex>);
  cDabij.contract(1.0, cVabij,"abij", cDabij,"abij", 0.0,"abij", fDivide);

  Scalar<complex> energy(*Cc4s::world);
  energy[""] = 2.0 * cDabij["abij"] * cVabij["abij"];
  complex direct(energy.get_val());
  LOG(1, "MP2") << "Direct energy computed analytically = " << direct << std::endl;

  energy[""] = -1.0 * cDabij["abij"] * cVabij["abji"];
  complex exchange(energy.get_val());
  LOG(1, "MP2") << "Exchange energy computed analytically = " << exchange << std::endl;
  return 2.0 * std::real(direct + exchange);
}

