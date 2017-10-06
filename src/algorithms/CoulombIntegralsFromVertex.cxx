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

  delete GammaGij;
  delete GammaGai;
  delete GammaGab;
  if (GammaGia) {
    delete GammaGia;
  }
}

void CoulombIntegralsFromVertex::dryRun() {
  // Read the Coulomb vertex GammaGqr
  DryTensor<complex> *GammaGqr(
    getTensorArgument<complex,DryTensor<complex>>("CoulombVertex")
  );

  // Read the Particle/Hole Eigenenergies
  DryTensor<> *epsi(
    getTensorArgument<double,DryTensor<double>>("HoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double,DryTensor<double>>("ParticleEigenEnergies")
  );

  LOG(0, "Integrals") <<
    "Reading Coulomb integrals form vertex " << GammaGqr->get_name() 
    << std::endl;

  // Compute the No,Nv,NG,Np
  int NG(GammaGqr->lens[0]);
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);

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

  bool realIntegrals = !getIntegerArgument("complex", 0);
  LOG(0, "CoulombIntegrals") << "Using "
    << (realIntegrals ? "real" : "complex") << " Coulomb integrals"
    << std::endl;

  // Allocate and compute GammaGab,GammaGai,GammaGij from GammaGqr
  int GaiLens[]   = {NG,Nv,No};
  int GiaLens[]   = {NG,No,Nv};
  int GabLens[]   = {NG,Nv,Nv};
  int GijLens[]   = {NG,No,No};

  DryTensor<complex> GammaGia(3, GiaLens, syms.data());
  DryTensor<complex> GammaGai(3, GaiLens, syms.data());
  DryTensor<complex> GammaGab(3, GabLens, syms.data());
  DryTensor<complex> GammaGij(3, GijLens, syms.data());

  if (realIntegrals) {
    dryCalculateRealIntegrals();
  } else {
    dryCalculateComplexIntegrals();
  }
}

void CoulombIntegralsFromVertex::calculateRealIntegrals() {
  Tensor<> *Vaibj(isArgumentGiven("PHPHCoulombIntegrals") ?
    new Tensor<>(4, vovo.data(), syms.data(), *Cc4s::world, "Vaibj") : nullptr);
  Tensor<> *Vabij(isArgumentGiven("PPHHCoulombIntegrals") ?
    new Tensor<>(4, vvoo.data(), syms.data(), *Cc4s::world, "Vabij") : nullptr);
  Tensor<> *Vijab(isArgumentGiven("HHPPCoulombIntegrals") ?
    new Tensor<>(4, oovv.data(), syms.data(), *Cc4s::world, "Vijab") : nullptr);
  Tensor<> *Vaijb(isArgumentGiven("PHHPCoulombIntegrals") ?
    new Tensor<>(4, voov.data(), syms.data(), *Cc4s::world, "Vaijb") : nullptr);
  Tensor<> *Vijkl(isArgumentGiven("HHHHCoulombIntegrals") ?
    new Tensor<>(4, oooo.data(), syms.data(), *Cc4s::world, "Vijkl") : nullptr);
  Tensor<> *Vijka(isArgumentGiven("HHHPCoulombIntegrals") ?
    new Tensor<>(4, ooov.data(), syms.data(), *Cc4s::world, "Vijka") : nullptr);
  Tensor<> *Vaijk(isArgumentGiven("PHHHCoulombIntegrals") ?
    new Tensor<>(4, vooo.data(), syms.data(), *Cc4s::world, "Vaijk") : nullptr);
  Tensor<> *Vabcd(isArgumentGiven("PPPPCoulombIntegrals") ?
    new Tensor<>(4, vvvv.data(), syms.data(), *Cc4s::world, "Vabcd") : nullptr);
  Tensor<> *Vabci(isArgumentGiven("PPPHCoulombIntegrals") ?
    new Tensor<>(4, vvvo.data(), syms.data(), *Cc4s::world, "Vabci") : nullptr);

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
  if (Vaibj) {
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vaibj->get_name() << std::endl;
    (*Vaibj)["aibj"]  = realGammaGab["Gab"] * realGammaGij["Gij"];
    (*Vaibj)["aibj"] += imagGammaGab["Gab"] * imagGammaGij["Gij"];
    allocatedTensorArgument("PHPHCoulombIntegrals", Vaibj);
  }
  if (Vaijb) {
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vaijb->get_name() << std::endl;
    (*Vaijb)["aijb"]  = realGammaGai["Gaj"] * realGammaGai["Gbi"];
    (*Vaijb)["aijb"] += imagGammaGai["Gaj"] * imagGammaGai["Gbi"];
    allocatedTensorArgument("PHHPCoulombIntegrals", Vaijb);
  }
  if (Vabij) {
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vabij->get_name() << std::endl;
    (*Vabij)["abij"]  = realGammaGai["Gai"] * realGammaGai["Gbj"];
    (*Vabij)["abij"] += imagGammaGai["Gai"] * imagGammaGai["Gbj"];
    allocatedTensorArgument("PPHHCoulombIntegrals", Vabij);
  }
  if (Vijab) {
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vijab->get_name() << std::endl;
    (*Vijab)["ijab"] = (*Vabij)["abij"];
    allocatedTensorArgument("HHPPCoulombIntegrals", Vijab);
  }
  if (Vijkl) {
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vijkl->get_name() << std::endl;
    (*Vijkl)["ijkl"]  = realGammaGij["Gik"] * realGammaGij["Gjl"];
    (*Vijkl)["ijkl"] += imagGammaGij["Gik"] * imagGammaGij["Gjl"];
    allocatedTensorArgument("HHHHCoulombIntegrals", Vijkl);
  }
  if (Vijka) {
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vijka->get_name() << std::endl;
    (*Vijka)["ijka"]  = realGammaGij["Gik"] * realGammaGai["Gaj"];
    (*Vijka)["ijka"] += imagGammaGij["Gik"] * imagGammaGai["Gaj"];
    allocatedTensorArgument("HHHPCoulombIntegrals", Vijka);
  }
  if (Vaijk) {
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vaijk->get_name() << std::endl;
    (*Vaijk)["aijk"]  = realGammaGai["Gaj"] * realGammaGij["Gik"];
    (*Vaijk)["aijk"] += imagGammaGai["Gaj"] * imagGammaGij["Gik"];
    allocatedTensorArgument("PHHHCoulombIntegrals", Vaijk);
  }
  if (Vabcd) {
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vabcd->get_name() << std::endl;
    (*Vabcd)["abcd"]  = realGammaGab["Gac"] * realGammaGab["Gbd"];
    (*Vabcd)["abcd"] += imagGammaGab["Gac"] * imagGammaGab["Gbd"];
    allocatedTensorArgument("PPPPCoulombIntegrals", Vabcd);
  }
  if (Vabci) {
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vabci->get_name() << std::endl;
    (*Vabci)["abci"]  = realGammaGab["Gac"] * realGammaGai["Gbi"];
    (*Vabci)["abci"] += imagGammaGab["Gac"] * imagGammaGai["Gbi"];
    allocatedTensorArgument("PPPHCoulombIntegrals", Vabci);
  }
}

void CoulombIntegralsFromVertex::dryCalculateRealIntegrals() {
  DryTensor<> *Vaibj(
    isArgumentGiven("PHPHCoulombIntegrals") ?
    new DryTensor<>(4, vovo.data(), syms.data()) : nullptr
  );
  DryTensor<> *Vabij(
    isArgumentGiven("PPHHCoulombIntegrals") ?
    new DryTensor<>(4, vvoo.data(), syms.data()) : nullptr
  );
  DryTensor<> *Vijab(
    isArgumentGiven("HHPPCoulombIntegrals") ?
    new DryTensor<>(4, oovv.data(), syms.data()) : nullptr
  );
  DryTensor<> *Vaijb(
    isArgumentGiven("PHHPCoulombIntegrals") ?
    new DryTensor<>(4, voov.data(), syms.data()) : nullptr
  );
  DryTensor<> *Vijkl(
    isArgumentGiven("HHHHCoulombIntegrals") ?
    new DryTensor<>(4, oooo.data(), syms.data()) : nullptr
  );
  DryTensor<> *Vijka(
    isArgumentGiven("HHHPCoulombIntegrals") ?
    new DryTensor<>(4, ooov.data(), syms.data()) : nullptr
  );
  DryTensor<> *Vaijk(
    isArgumentGiven("PHHHCoulombIntegrals") ?
    new DryTensor<>(4, vooo.data(), syms.data()) : nullptr
  );
  DryTensor<> *Vabcd(
    isArgumentGiven("PPPPCoulombIntegrals") ?
    new DryTensor<>(4, vvvv.data(), syms.data()) : nullptr
  );
  DryTensor<> *Vabci(
    isArgumentGiven("PPPHCoulombIntegrals") ?
    new DryTensor<>(4, vvvo.data(), syms.data()) : nullptr
  );

  if (Vaibj) {
    allocatedTensorArgument("PHPHCoulombIntegrals", Vaibj);
  }
  if (Vabij) {
    allocatedTensorArgument("PPHHCoulombIntegrals", Vabij);
  }
  if (Vijab) {
    allocatedTensorArgument("HHPPCoulombIntegrals", Vijab);
  }
  if (Vaijb) {
    allocatedTensorArgument("PHHPCoulombIntegrals", Vaijb);
  }
  if (Vijkl) {
    allocatedTensorArgument("HHHHCoulombIntegrals", Vijkl);
  }
  if (Vijka) {
    allocatedTensorArgument("HHHPCoulombIntegrals", Vijka);
  }
  if (Vaijk) {
    allocatedTensorArgument("PHHHCoulombIntegrals", Vaijk);
  }
  if (Vabcd) {
    allocatedTensorArgument("PPPPCoulombIntegrals", Vabcd);
  }
  if (Vabci) {
    allocatedTensorArgument("PPPHCoulombIntegrals", Vabci);
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
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vabij->get_name() << std::endl;
    (*Vabij)["abij"] = conjTransposeGammaGai["Gai"] * (*GammaGai)["Gbj"];
    allocatedTensorArgument<complex>("PPHHCoulombIntegrals", Vabij);
  }

  if (Vaijb) {
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vaijb->get_name() << std::endl;
    (*Vaijb)["aijb"] = conjTransposeGammaGai["Gaj"] * (*GammaGia)["Gib"];
    allocatedTensorArgument<complex>("PHHPCoulombIntegrals", Vaijb);
  }

  if (Vijab) {
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vijab->get_name() << std::endl;
    (*Vijab)["ijab"] = conjTransposeGammaGia["Gia"] * (*GammaGia)["Gjb"];
    allocatedTensorArgument<complex>("HHPPCoulombIntegrals", Vijab);
  }

  if (Vaibj) {
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vaibj->get_name() << std::endl;
    (*Vaibj)["aibj"] = conjTransposeGammaGab["Gab"] * (*GammaGij)["Gij"];
    allocatedTensorArgument<complex>("PHPHCoulombIntegrals", Vaibj);
  }

  if (Vijkl) {
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vijkl->get_name() << std::endl;
    (*Vijkl)["ijkl"] = conjTransposeGammaGij["Gik"] * (*GammaGij)["Gjl"];
    allocatedTensorArgument<complex>("HHHHCoulombIntegrals", Vijkl);
  }

  if (Vijka) {
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vijka->get_name() << std::endl;
    (*Vijka)["ijka"] = conjTransposeGammaGij["Gik"] * (*GammaGia)["Gja"];
    allocatedTensorArgument<complex>("HHHPCoulombIntegrals", Vijka);
  }

  if (Vaijk) {
    LOG(1, "CoulombIntegrals") << "Evaluating "
                               << Vaijk->get_name() << std::endl;
    (*Vaijk)["aijk"] = conjTransposeGammaGai["Gaj"] * (*GammaGij)["Gik"];
    allocatedTensorArgument<complex>("PHHHCoulombIntegrals", Vaijk);
  }
}

void CoulombIntegralsFromVertex::dryCalculateComplexIntegrals() {
  DryTensor<complex> *Vaibj(
    isArgumentGiven("PHPHCoulombIntegrals") ?
    new DryTensor<complex>(4, vovo.data(), syms.data())
    : nullptr
  );
  DryTensor<complex> *Vabij(
    isArgumentGiven("PPHHCoulombIntegrals") ?
    new DryTensor<complex>(4, vvoo.data(), syms.data())
    : nullptr
  );
  DryTensor<complex> *Vijab(
    isArgumentGiven("HHPPCoulombIntegrals") ?
    new DryTensor<complex>(4, oovv.data(), syms.data())
    : nullptr
  );
  DryTensor<complex> *Vaijb(
    isArgumentGiven("PHHPCoulombIntegrals") ?
    new DryTensor<complex>(4, voov.data(), syms.data())
    : nullptr
  );
  DryTensor<complex> *Vijkl(
    isArgumentGiven("HHHHCoulombIntegrals") ?
    new DryTensor<complex>(4, oooo.data(), syms.data())
    : nullptr
  );
  DryTensor<complex> *Vijka(
    isArgumentGiven("HHHPCoulombIntegrals") ?
    new DryTensor<complex>(4, ooov.data(), syms.data())
    : nullptr
  );
  DryTensor<complex> *Vaijk(
    isArgumentGiven("PHHHCoulombIntegrals") ?
    new DryTensor<complex>(4, vooo.data(), syms.data())
    : nullptr
  );

  if (Vaibj) {
    allocatedTensorArgument<complex,DryTensor<complex>>("PHPHCoulombIntegrals", Vaibj);
  }
  if (Vabij) {
    allocatedTensorArgument<complex,DryTensor<complex>>("PPHHCoulombIntegrals", Vabij);
  }
  if (Vijab) {
    allocatedTensorArgument<complex,DryTensor<complex>>("HHPPCoulombIntegrals", Vijab);
  }
  if (Vaijb) {
    allocatedTensorArgument<complex,DryTensor<complex>>("PHHPCoulombIntegrals", Vaijb);
  }
  if (Vijkl) {
    allocatedTensorArgument<complex,DryTensor<complex>>("HHHHCoulombIntegrals", Vijkl);
  }
  if (Vijka) {
    allocatedTensorArgument<complex,DryTensor<complex>>("HHHPCoulombIntegrals", Vijka);
  }
  if (Vaijk) {
    allocatedTensorArgument<complex,DryTensor<complex>>("PHHHCoulombIntegrals", Vaijk);
  }
}

