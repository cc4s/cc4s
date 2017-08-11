#include <algorithms/CoulombIntegralsFromVertex.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(CoulombIntegralsFromVertex);

CoulombIntegralsFromVertex::CoulombIntegralsFromVertex(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

CoulombIntegralsFromVertex::~CoulombIntegralsFromVertex() {
}

void CoulombIntegralsFromVertex::run() {
  // Read the Coulomb vertex GammaGqr
  Tensor<complex> *GammaGqr( getTensorArgument<complex>("CoulombVertex"));

  // Read the Particle/Hole Eigenenergies
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));

  LOG(0, "Integrals") <<
    "Reading Coulomb integrals form vertex " << GammaGqr->get_name() 
    << std::endl;

  // Compute the No,Nv,NG,Np
  int NG(GammaGqr->lens[0]);
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np(GammaGqr->lens[1]);

  // Allocate coulomb integrals Vabij Vaibj Vaijb Vijkl Vabcd
  syms = std::array<int,4>{{ NS, NS, NS, NS }};
  vvvv = std::array<int,4>{{ Nv, Nv, Nv, Nv }};
  vovo = std::array<int,4>{{ Nv, No, Nv, No }};
  vvoo = std::array<int,4>{{ Nv, Nv, No, No }};
  voov = std::array<int,4>{{ Nv, No, No, Nv }};
  oovv = std::array<int,4>{{ No, No, Nv, Nv }};
  oooo = std::array<int,4>{{ No, No, No, No }};
  ooov = std::array<int,4>{{ No, No, No, Nv }};
  vooo = std::array<int,4>{{ Nv, No, No, No }};
  vvvo = std::array<int,4>{{ Nv, Nv, Nv, No }};

  /*
  Tensor<complex> diffGammaGqr(*GammaGqr);
  diffGammaGqr["Gqr"] -= (*GammaGqr)["Grq"];
  double diffGamma(frobeniusNorm(diffGammaGqr));
  LOG(1, "CoulombIntegrals") << "|GammaGqr-GammaGrq|=" << diffGamma
  << std::endl;
  */

  // diffGamma = 0 iff orbitals are real valued
  //  double const threshold(1e-10);
  //  bool realIntegrals(diffGamma < threshold);

  bool realIntegrals = !getIntegerArgument("complex", 0);
  LOG(0, "CoulombIntegrals") << "Using "
    << (realIntegrals ? "real" : "complex") << " Coulomb integrals"
    << std::endl;

  int aStart(Np-Nv), aEnd(Np);
  int iStart(0), iEnd(No);
  int GijStart[] = {0, iStart,iStart};
  int GijEnd[]   = {NG,iEnd,  iEnd};
  int GiaStart[] = {0, iStart,aStart};
  int GiaEnd[]   = {NG,iEnd,  aEnd};
  int GaiStart[] = {0, aStart,iStart};
  int GaiEnd[]   = {NG,aEnd,  iEnd};
  int GabStart[] = {0, aStart,aStart};
  int GabEnd[]   = {NG,aEnd,  aEnd};
  GammaGij = new Tensor<complex>(GammaGqr->slice(GijStart, GijEnd));
  GammaGia = realIntegrals ?
    nullptr : new Tensor<complex>(GammaGqr->slice(GiaStart, GiaEnd));
  GammaGai = new Tensor<complex>(GammaGqr->slice(GaiStart, GaiEnd));
  GammaGab = new Tensor<complex>(GammaGqr->slice(GabStart, GabEnd));
  

  if (realIntegrals) {
    calculateRealIntegrals();
  } else {
    calculateComplexIntegrals();
  }
}

// FIXME: update dryRun to work in the complex case as well
void CoulombIntegralsFromVertex::dryRun() {
  // Read the Coulomb vertex GammaGqr
  DryTensor<complex> *GammaGqr(getTensorArgument<complex, 
                               DryTensor<complex>>("CoulombVertex"));

  // Read the Particle/Hole Eigenenergies
  DryTensor<> *epsi(getTensorArgument
                    <double, DryTensor<double>>("HoleEigenEnergies"));
  DryTensor<> *epsa(getTensorArgument
                    <double, DryTensor<double>>("ParticleEigenEnergies"));

  // Compute the No,Nv,NG
  int NG(GammaGqr->lens[0]);
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

  // Allocate coulomb integrals Vabij Vaibj Vaijb Vijkl Vabcd
  int syms[] = { NS, NS, NS, NS };
  int vvvv[] = { Nv, Nv, Nv, Nv };
  int vovo[] = { Nv, No, Nv, No };
  int vvoo[] = { Nv, Nv, No, No };
  int oooo[] = { No, No, No, No };
  int ooov[] = { No, No, No, Nv };
  int vvvo[] = { Nv, Nv, Nv, No };

  DryTensor<> *Vabcd(isArgumentGiven("PPPPCoulombIntegrals") 
                     ?new DryTensor<>(4, vvvv, syms) : nullptr);
  DryTensor<> *Vaibj(isArgumentGiven("PHPHCoulombIntegrals") 
                     ?new DryTensor<>(4, vovo, syms) : nullptr);
  DryTensor<> *Vabij(isArgumentGiven("PPHHCoulombIntegrals") ?
                     new DryTensor<>(4, vvoo, syms) : nullptr);
  DryTensor<> *Vijkl(isArgumentGiven("HHHHCoulombIntegrals") ?
                     new DryTensor<>(4, oooo, syms) : nullptr);
  DryTensor<> *Vijka(isArgumentGiven("HHHPCoulombIntegrals") ?
                     new DryTensor<>(4, ooov, syms) : nullptr);
  DryTensor<> *Vabci(isArgumentGiven("PPPHCoulombIntegrals") ?
                     new DryTensor<>(4, vvvo, syms) : nullptr);

  if (Vabcd) {
    allocatedTensorArgument("PPPPCoulombIntegrals", Vabcd);
  }
  if (Vaibj) {
    allocatedTensorArgument("PHPHCoulombIntegrals", Vaibj);
  }
  if (Vabij) {
    allocatedTensorArgument("PPHHCoulombIntegrals", Vabij);
  }
  if (Vijkl) {
    allocatedTensorArgument("HHHHCoulombIntegrals", Vijkl);
  }
  if (Vijka) {
    allocatedTensorArgument("HHHPCoulombIntegrals", Vijka);
  }
  if (Vabci) {
    allocatedTensorArgument("PPPHCoulombIntegrals", Vabci);
  }

  // Allocate and compute GammaGab,GammaGai,GammaGij from GammaGqr
  int GaiLens[]   = {NG,Nv,No};
  int GabLens[]   = {NG,Nv,Nv};
  int GijLens[]   = {NG,No,No};

  DryTensor<complex> GammaGai(3, GaiLens, syms);
  DryTensor<complex> GammaGab(3, GabLens, syms);
  DryTensor<complex> GammaGij(3, GijLens, syms);

  // Split GammaGab,GammaGai,GammaGij into real and imaginary parts
  DryTensor<> realGammaGai(3, GaiLens, syms);
  DryTensor<> imagGammaGai(3, GaiLens, syms);

  DryTensor<> realGammaGab(3, GabLens, syms);
  DryTensor<> imagGammaGab(3, GabLens, syms);

  DryTensor<> realGammaGij(3, GijLens, syms);
  DryTensor<> imagGammaGij(3, GijLens, syms);
}

void CoulombIntegralsFromVertex::calculateRealIntegrals() {
  Tensor<> *Vabcd(isArgumentGiven("PPPPCoulombIntegrals") ?
    new Tensor<>(4, vvvv.data(), syms.data(), *Cc4s::world, "Vabcd") : nullptr);
  Tensor<> *Vaibj(isArgumentGiven("PHPHCoulombIntegrals") ?
    new Tensor<>(4, vovo.data(), syms.data(), *Cc4s::world, "Vaibj") : nullptr);
  Tensor<> *Vabij(isArgumentGiven("PPHHCoulombIntegrals") ?
    new Tensor<>(4, vvoo.data(), syms.data(), *Cc4s::world, "Vabij") : nullptr);
  Tensor<> *Vijkl(isArgumentGiven("HHHHCoulombIntegrals") ?
    new Tensor<>(4, oooo.data(), syms.data(), *Cc4s::world, "Vijkl") : nullptr);
  Tensor<> *Vijka(isArgumentGiven("HHHPCoulombIntegrals") ?
    new Tensor<>(4, ooov.data(), syms.data(), *Cc4s::world, "Vijka") : nullptr);
  Tensor<> *Vabci(isArgumentGiven("PPPHCoulombIntegrals") ?
    new Tensor<>(4, vvvo.data(), syms.data(), *Cc4s::world, "Vabci") : nullptr);

  if (Vabcd) {
    allocatedTensorArgument("PPPPCoulombIntegrals", Vabcd);
  }
  if (Vaibj) {
    allocatedTensorArgument("PHPHCoulombIntegrals", Vaibj);
  }
  if (Vabij) {
    allocatedTensorArgument("PPHHCoulombIntegrals", Vabij);
  }
  if (Vijkl) {
    allocatedTensorArgument("HHHHCoulombIntegrals", Vijkl);
  }
  if (Vijka) {
    allocatedTensorArgument("HHHPCoulombIntegrals", Vijka);
  }
  if (Vabci) {
    allocatedTensorArgument("PPPHCoulombIntegrals", Vabci);
  }

  // Split GammaGab,GammaGai,GammaGia,GammaGij into real and imaginary parts
  Tensor<> realGammaGai(3, GammaGai->lens, GammaGai->sym, 
                        *GammaGai->wrld, "RealGammaGai");
  Tensor<> imagGammaGai(3, GammaGai->lens, GammaGai->sym, 
                        *GammaGai->wrld, "ImagGammaGai");
  fromComplexTensor(*GammaGai, realGammaGai, imagGammaGai);

  Tensor<> realGammaGab(3, GammaGab->lens, GammaGab->sym, 
                        *GammaGab->wrld, "RealGammaGab");
  Tensor<> imagGammaGab(3, GammaGab->lens, GammaGab->sym, 
                        *GammaGab->wrld, "ImagGammaGab");
  fromComplexTensor(*GammaGab, realGammaGab, imagGammaGab);

  Tensor<> realGammaGij(3, GammaGij->lens, GammaGij->sym, 
                        *GammaGij->wrld, "RealGammaGij");
  Tensor<> imagGammaGij(3, GammaGij->lens, GammaGij->sym, 
                        *GammaGij->wrld, "ImagGammaGij");
  fromComplexTensor(*GammaGij, realGammaGij, imagGammaGij);

  // Compute the integrals Vabij Vaibj Vaijb Vijkl Vabcd
  if (Vabcd) {
    LOG(1, "CoulombIntegrals") << "Evaluating " 
                               << Vabcd->get_name() << std::endl;
    (*Vabcd)["abcd"]  = realGammaGab["Gac"] * realGammaGab["Gbd"];
    (*Vabcd)["abcd"] += imagGammaGab["Gac"] * imagGammaGab["Gbd"];
  }
  if (Vaibj) {
    LOG(1, "CoulombIntegrals") << "Evaluating " 
                               << Vaibj->get_name() << std::endl;
    (*Vaibj)["aibj"]  = realGammaGab["Gab"] * realGammaGij["Gij"];
    (*Vaibj)["aibj"] += imagGammaGab["Gab"] * imagGammaGij["Gij"];
  }
  if (Vabij) {
    LOG(1, "CoulombIntegrals") << "Evaluating " 
                               << Vabij->get_name() << std::endl;
    (*Vabij)["abij"]  = realGammaGai["Gai"] * realGammaGai["Gbj"];
    (*Vabij)["abij"] += imagGammaGai["Gai"] * imagGammaGai["Gbj"];
    double r(frobeniusNorm(*Vabij));
    LOG(1, "CoulombIntegrals") << "|Vabij|=" << r << std::endl;
  }
  if (Vijkl) {
    LOG(1, "CoulombIntegrals") << "Evaluating " 
                               << Vijkl->get_name() << std::endl;
    (*Vijkl)["ijkl"]  = realGammaGij["Gik"] * realGammaGij["Gjl"];
    (*Vijkl)["ijkl"] += imagGammaGij["Gik"] * imagGammaGij["Gjl"];
  }
  if (Vijka) {
    LOG(1, "CoulombIntegrals") << "Evaluating " 
                               << Vijka->get_name() << std::endl;
    (*Vijka)["ijka"]  = realGammaGij["Gik"] * realGammaGai["Gaj"];
    (*Vijka)["ijka"] += imagGammaGij["Gik"] * imagGammaGai["Gaj"];
  }
  if (Vabci) {
    LOG(1, "CoulombIntegrals") << "Evaluating " 
                               << Vabci->get_name() << std::endl;
    (*Vabci)["abci"]  = realGammaGab["Gac"] * realGammaGai["Gbi"];
    (*Vabci)["abci"] += imagGammaGab["Gac"] * imagGammaGai["Gbi"];
  }
}

void CoulombIntegralsFromVertex::calculateComplexIntegrals() {
  Tensor<complex> *Vabij(
    isArgumentGiven("PPHHCoulombIntegrals") ?
      new Tensor<complex>(4, vvoo.data(), syms.data(), *Cc4s::world, "Vabij") :
      nullptr
  );

  Tensor<complex> *Vijab(
    // TODO: HHPP is always conj(PPHH)
    isArgumentGiven("HHPPCoulombIntegrals") ?
      new Tensor<complex>(4, oovv.data(), syms.data(), *Cc4s::world, "Vijab") :
      nullptr
  );

  Tensor<complex> *Vaijb(
    isArgumentGiven("PHHPCoulombIntegrals") ?
      new Tensor<complex>(4, voov.data(), syms.data(), *Cc4s::world, "Vaijb") :
      nullptr
  );

  
  Tensor<complex> *Vaibj(
    isArgumentGiven("PHPHCoulombIntegrals") ?
      new Tensor<complex>(4, vovo.data(), syms.data(), *Cc4s::world, "Vaibj") :
      nullptr
  );

  Tensor<complex> *Vijkl(
    isArgumentGiven("HHHHCoulombIntegrals") ?
      new Tensor<complex>(4, oooo.data(), syms.data(), *Cc4s::world, "Vijkl") :
      nullptr
  );
  
  Tensor<complex> *Vijka(
    isArgumentGiven("HHHPCoulombIntegrals") ?
      new Tensor<complex>(4, ooov.data(), syms.data(), *Cc4s::world, "Vijka") :
      nullptr
  );

  Tensor<complex> *Vaijk(
    // TODO: PHHH is always conj(Permute(HHHP))
    isArgumentGiven("PHHHCoulombIntegrals") ?
      new Tensor<complex>(4, vooo.data(), syms.data(), *Cc4s::world, "Vaijk") :
      nullptr
  );

  Univar_Function<complex> fConj(conj<complex>);

  Tensor<complex> conjTransposeGammaGai(false, *GammaGai);
  conjTransposeGammaGai.sum(1.0,*GammaGia,"Gia", 0.0,"Gai", fConj);

  Tensor<complex> conjTransposeGammaGia(false, *GammaGia);
  conjTransposeGammaGia.sum(1.0,*GammaGai,"Gai", 0.0,"Gia", fConj);

  Tensor<complex> conjTransposeGammaGij(false, *GammaGij);
  conjTransposeGammaGij.sum(1.0,*GammaGij,"Gji", 0.0,"Gij", fConj);

  Tensor<complex> conjTransposeGammaGab(false, *GammaGab);
  conjTransposeGammaGab.sum(1.0,*GammaGab,"Gba", 0.0,"Gab", fConj);


  if (Vabij) {
    (*Vabij)["abij"] = conjTransposeGammaGai["Gai"] * (*GammaGai)["Gbj"];
  }

  if (Vaijb) {
    (*Vaijb)["aijb"] = conjTransposeGammaGai["Gaj"] * (*GammaGia)["Gib"];
  }

  if (Vijab) {
    (*Vijab)["ijab"] = conjTransposeGammaGia["Gia"] * (*GammaGia)["Gjb"];
  }

  if (Vaibj) {
    (*Vaibj)["aibj"] = conjTransposeGammaGab["Gab"] * (*GammaGij)["Gij"];
  }

  if (Vijkl) {
    (*Vijkl)["ijkl"] = conjTransposeGammaGij["Gik"] * (*GammaGij)["Gjl"];
  }

  if (Vijka) {
    (*Vijka)["ijka"] = conjTransposeGammaGij["Gik"] * (*GammaGia)["Gja"];
  }

  if (Vaijk) {
    (*Vaijk)["aijk"] = conjTransposeGammaGai["Gaj"] * (*GammaGij)["Gik"];
  }

  /*
  Tensor<> realVabij(4, Vabij->lens, Vabij->sym, *Vabij->wrld, "realVabij");
  Tensor<> imagVabij(4, Vabij->lens, Vabij->sym, *Vabij->wrld, "imagVabij");
  fromComplexTensor(*Vabij, realVabij, imagVabij);
  double r(frobeniusNorm(realVabij));
  double i(frobeniusNorm(imagVabij));
  LOG(1, "CoulombIntegrals") << "|Re(Vabij)|=" << r << ", |Im(Vabij)|=" << i << std::endl;
  */

  if (Vabij) {
    allocatedTensorArgument<complex>("PPHHCoulombIntegrals", Vabij);
  }
  if (Vaijb) {
    allocatedTensorArgument<complex>("PHHPCoulombIntegrals", Vaijb);
  }
  if (Vijab) {
    allocatedTensorArgument<complex>("HHPPCoulombIntegrals", Vijab);
  }
  if (Vaibj) {
    allocatedTensorArgument<complex>("PHPHCoulombIntegrals", Vaibj);
  }
  if (Vijkl) {
    allocatedTensorArgument<complex>("HHHHCoulombIntegrals", Vijkl);
  }
  if (Vijka) {
    allocatedTensorArgument<complex>("HHHPCoulombIntegrals", Vijka);
  }
  if (Vaijk) {
    allocatedTensorArgument<complex>("PHHHCoulombIntegrals", Vaijk);
  }
}

