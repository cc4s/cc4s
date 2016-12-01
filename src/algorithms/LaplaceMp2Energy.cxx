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
}

void LaplaceMp2Energy::run() {
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<> *Tn  (getTensorArgument("ImaginaryTimePoints"));
  Tensor<> *Wn  (getTensorArgument("ImaginaryTimeWeights"));

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

  Hin["in"]  = (*epsi)["i"];
  Hin["in"] += (-1.0) * mu;
  Hin["in"] *= (*Tn)["n"];

  Hin["in"] = CTF::Function<double, double>(
    std::function<double(double)>([](double c){ return std::exp(c); })
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
  int NG(LambdaGR->lens[0]);
  int  RR[] =     { NR, NR };
  int RRn[] = { NR, NR, Nn };

  Tensor<complex> VRS  (2, RR,  syms, *PirR->wrld,   "VRS");
  Tensor<complex> GpRSn(2, RRn, syms, *PirR->wrld, "GpRSn");
  Tensor<complex> GhRSn(2, RRn, syms, *PirR->wrld, "GhRSn");

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
  VRS["RS"] = conjLambdaGR["GR"] * (*LambdaGR)["GS"];

  Tensor<complex> conjPiiR(false, PiiR);
  conjPiiR.set_name("ConjPiiR");
  conjPiiR.sum(1.0, PiiR,"iR", 0.0,"iR", fConj);
  GhRSn["RSn"] = conjPiiR["iS"] * PiiR["iR"] * Hin["in"];

  Tensor<complex> conjPiaR(false, PiaR);
  conjPiaR.set_name("ConjPiaR");
  conjPiaR.sum(1.0, PiaR,"aR", 0.0,"aR", fConj);
  GpRSn["RSn"] = conjPiaR["aS"] * PiaR["aR"] * Pan["an"];

  Scalar<> energy(*Cc4s::world);
  double e, dire, exce;

  energy[""] = VRS["RS"] * GpRSn["RTn"] * GhRSn["TRn"] * GpRSn["SUn"] * GhRSn["USn"] * VRS["TU"];
  dire =  2.0 * energy.get_val();
  energy[""] =  VRS["RS"] * GpRSn["SUn"] * GhRSn["TSn"] * VRS["TU"] * GpRSn["RTn"] * GhRSn["URn"];
  exce = -1.0 * energy.get_val();
  e = dire + exce;

  LOG(0, "MP2") << "e=" << e << std::endl;
  LOG(1, "MP2") << "MP2d=" << dire << std::endl;
  LOG(1, "MP2") << "MP2x=" << exce << std::endl;

  setRealArgument("Mp2Energy", e);
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

  int Np(PirR->lens[0]);
  int NR(PirR->lens[1]);
  int NG(LambdaGR->lens[0]);
  int RR[] =   { NR, NR };
  int RRn[] =  { NR, NR, Nn };

  DryTensor<complex> VRS  (2, RR,  syms, SOURCE_LOCATION);
  DryTensor<complex> GpRSn(2, RRn, syms, SOURCE_LOCATION);
  DryTensor<complex> GhRSn(2, RRn, syms, SOURCE_LOCATION);

  // Allocate and compute PiaR
  int vR[] =   { Nv, NR };
  DryTensor<complex> PiaR(2, vR,  syms, SOURCE_LOCATION);

  // Allocate and compute PiiR
  int oR[] =   { No, NR };
  DryTensor<complex> PiiR(2, oR,  syms, SOURCE_LOCATION);

  DryTensor<complex> conjLambdaGR(*LambdaGR);

  DryTensor<complex> conjPiiR(PiiR);

  DryTensor<complex> conjPiaR(PiaR);

  DryScalar<> energy();
}

