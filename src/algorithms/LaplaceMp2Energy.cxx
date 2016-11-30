#include <algorithms/LaplceMp2Energy.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(LaplceMp2Energy);

LaplceMp2Energy::LaplceMp2Energy(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

LaplceMp2Energy::~LaplceMp2Energy() {
}

void LaplceMp2Energy::run() {
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  Tensor<> *Tn  (getTensorArgument("LaplaceGridPoints"));
  Tensor<> *Wn  (getTensorArgument("LaplaceWeights"));

  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Nn(Tn->lens[0]);
  int on[] = { No, Nn };
  int vn[] = { Nv, Nn };
  Tensor<> Pan(2, vn, syms, *epsi->wrld, "Pan");
  Tensor<> Hin(2, on, syms, *epsi->wrld, "Hin");

  Tensor<complex> *PirR(getTensorArgument<complex>("FactorOrbitals"));
  PirR->set_name("PirR");
  Tensor<complex> *LambdaGR(getTensorArgument<complex>("CoulombFactors"));
  LambdaGR->set_name("LambdaGR"); 

  int Np(PirR->lens[0]);
  int NR(PirR->lens[1]);
  int NG(LambdaGR->lens[0]);
  int RR[] =   { Rx, Ry };
  int RRn[] =  { Rx, Ry, Nn };
  int syms[] = { NS, NS, NS, NS };

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
  Tensor<complex> conjLeftLambdaGR(false, leftLambdaGR);
  conjLeftLambdaGR.set_name("ConjLeftLambdaGR");
  conjLeftLambdaGR.sum(1.0, leftLambdaGR,"GR", 0.0,"GR", fConj);
  VRS["RS"] = conjLeftLambdaGR["GR"] * rightLambdaGR["GS"];

  Univar_Function<complex> fConj(&cc4s::conj<complex>);
  Tensor<complex> conjPiiR(false, PiiR);
  conjPiiR.set_name("ConjPiiR");
  conjPiiR.sum(1.0, PiiR,"iR", 0.0,"iR", fConj);
  GhRSn["RSn"] = conjPiiR["iS"] * PiiR["iR"] * Hin["in"];

  Univar_Function<complex> fConj(&cc4s::conj<complex>);
  Tensor<complex> conjPiaR(false, PiaR);
  conjPiaR.set_name("ConjPiaR");
  conjPiaR.sum(1.0, PiaR,"aR", 0.0,"aR", fConj);
  GpRSn["RSn"] = conjPiaR["aS"] * PiaR["aR"] * Pan["an"];

  Scalar<> energy(*Cc4s::world);
  double e, dire, exce;

  energy[""] = GpRSn["RTn"] * GhRSn["TRn"] * VRS["RS"] * GpRSn["SUn"] * GhRSn["USn"] * VRS["TU"];
  dire =  2.0 * energy.get_val();
  energy[""] = GpRSn["RTn"] * GhRSn["URn"] * VRS["RS"] * GpRSn["SUn"] * GhRSn["RSn"] * VRS["TU"];
  exce = -1.0 * energy.get_val();
  e = dire + exce;

  LOG(0, "MP2") << "e=" << e << std::endl;
  LOG(1, "MP2") << "MP2d=" << dire << std::endl;
  LOG(1, "MP2") << "MP2x=" << exce << std::endl;

  setRealArgument("Mp2Energy", e);
}

void LaplceMp2Energy::dryRun() {
  //DryTensor<> *Vabij(
  getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals");
  //);

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

  // Allocate the doubles amplitudes
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { Nv, Nv, No, No };
  DryTensor<> Tabij(4, vvoo, syms);

  DryScalar<> energy();
}

