#include <algorithms/LaplaceMp2Energy.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <Options.hpp>

using namespace CTF;
using namespace cc4s;

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
  Tensor<> *Tn  (getTensorArgument("ImaginaryTimePoints"));

  // TODO: slice over imaginary time samples to conserve memory
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Nn(Tn->lens[0]);
  int on[] = { No, Nn };
  int vn[] = { Nv, Nn };
  int syms[] = { NS, NS, NS, NS };
  Tensor<> Pan(2, vn, syms, *epsi->wrld, "Pan");
  Tensor<> Hin(2, on, syms, *epsi->wrld, "Hin");

  double mu(getRealArgument("ChemicalPotential"));
  Pan["an"]  = (*epsa)["a"];
  Pan["an"] += (-1.0) * mu;
  Pan["an"] *= (*Tn)["n"];

  Pan["an"] = CTF::Function<double, double>(
    std::function<double(double)>([](double c){ return std::exp(-c); })
  ) (
    Pan["an"]
  );
  LOG(0, "IMP2") << "Particle time propagator set up" << std::endl;

  Hin["in"]  = (*epsi)["i"];
  Hin["in"] += (-1.0) * mu;
  Hin["in"] *= (*Tn)["n"];

  Hin["in"] = CTF::Function<double, double>(
    std::function<double(double)>([](double c){ return std::exp(c); })
  ) (
    Hin["in"]
  );
  LOG(0, "IMP2") << "Hole time propagator set up" << std::endl;

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

  VRS =   new Tensor<complex>(2, RR,  syms, *PirR->wrld,   "VRS");
  GpRSn = new Tensor<complex>(3, RRn, syms, *PirR->wrld, "GpRSn");
  GhRSn = new Tensor<complex>(3, RRn, syms, *PirR->wrld, "GhRSn");

  // Allocate and compute PiaR
  int aRStart[] = {No , 0};
  int aREnd[]   = {Np ,NR};
  Tensor<complex> PiaR(PirR->slice(aRStart,aREnd));
  PiaR.set_name("PiaR");

  // Allocate and compute PiiR
  int iRStart[] = {0 , 0};
  int iREnd[]   = {No ,NR};
  Tensor<complex> PiiR(PirR->slice(iRStart,iREnd));
  PiiR.set_name("PiiR");

  Univar_Function<complex> fConj(&cc4s::conj<complex>);
  Tensor<complex> conjLambdaGR(false, LambdaGR);
  conjLambdaGR.set_name("ConjLambdaGR");
  conjLambdaGR.sum(1.0, *LambdaGR,"GR", 0.0,"GR", fConj);
  (*VRS)["RS"] = conjLambdaGR["GR"] * (*LambdaGR)["GS"];
  // FIXME: imaginary part manually discarded, should be zero
  VRS->sum(-0.5, *VRS,"RS", 0.5,"RS", fConj);
  LOG(0, "IMP2") << "Coulomb propagator set up" << std::endl;

  // TODO: check for real factor orbitals
  (*GhRSn)["RSn"] = PiiR["iS"] * PiiR["iR"] * Hin["in"];
  LOG(0, "IMP2") << "Hole propagator set up" << std::endl;

  (*GpRSn)["RSn"] = PiaR["aS"] * PiaR["aR"] * Pan["an"];
  LOG(0, "IMP2") << "Particle propagator set up" << std::endl;


  // get numerical weights
  Tensor<> *realwn(getTensorArgument("ImaginaryTimeWeights"));
  wn = new Tensor<complex>(1, &Nn, syms, *PirR->wrld, "wn");
  toComplexTensor(*realwn, *wn);


  double dire(calculateDirectTerm());
//  double exce(calculateExchangeTerm());
//  double e(dire + exce);

//  LOG(0, "MP2") << "e=" << e << std::endl;
  LOG(1, "MP2") << "MP2d=" << dire << std::endl;
//  LOG(1, "MP2") << "MP2x=" << exce << std::endl;

  setRealArgument("Mp2Energy", dire);
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
  ChiVRSn["RSn"] *= (*VRS)["RS"];
  energy[""] = (*wn)["n"] * ChiVRSn["RSn"] * ChiVRSn["SRn"];
  LOG(0, "IMP2") << "Direct energy computed" << std::endl;
  return 2.0 * std::real(energy.get_val());
}

double LaplaceMp2Energy::calculateExchangeTerm() {
//  energy[""] =  VRS["RS"] * GpRSn["SUn"] * GhRSn["TSn"] * VRS["TU"] * GpRSn["RTn"] * GhRSn["URn"] * (*Wn)["n"];
//  LOG(0, "IMP2") << "Exchange energy computed" << std::endl;
//  exce = -1.0 * energy.get_val();
  return 0.0;
}

