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
  // hole and particle states may overlap at finite temperature
  int Np(GammaGqr->lens[1]);
  int iStart(0), iEnd(No);
  int aStart(Np-Nv), aEnd(Np);

  // Allocate coulomb integrals Vabij Vaibj Vaijb Vijkl Vabcd
  int syms[] = { NS, NS, NS, NS };
  int vvvv[] = { Nv, Nv, Nv, Nv };
  int vovo[] = { Nv, No, Nv, No };
  int vvoo[] = { Nv, Nv, No, No };
  int oooo[] = { No, No, No, No };
  int ooov[] = { No, No, No, Nv };
  int vvvo[] = { Nv, Nv, Nv, No };

  // Indices for integrals created from already existing
  int oovv[] = { No, No, Nv, Nv };
  int ovoo[] = { No, Nv, No, No };
  int ovov[] = { No, Nv, No, Nv };
  int voov[] = { Nv, No, No, Nv };
  int ovvv[] = { No, Nv, Nv, Nv };
  int vvov[] = { Nv, Nv, No, Nv };

  int antisymmetrize(getIntegerArgument("antisymmetrize", 0));

  Tensor<> *Vabcd(
    isArgumentGiven("PPPPCoulombIntegrals") ? 
    new Tensor<>(4, vvvv, syms, *Cc4s::world, "Vabcd") : nullptr
  );
  Tensor<> *Vaibj(
    isArgumentGiven("PHPHCoulombIntegrals") ?
    new Tensor<>(4, vovo, syms, *Cc4s::world, "Vaibj") : nullptr
  );
  Tensor<> *Vabij(
    isArgumentGiven("PPHHCoulombIntegrals") ?
    new Tensor<>(4, vvoo, syms, *Cc4s::world, "Vabij") : nullptr
  );
  Tensor<> *Vijkl(
    isArgumentGiven("HHHHCoulombIntegrals") ?
    new Tensor<>(4, oooo, syms, *Cc4s::world, "Vijkl") : nullptr
  );
  Tensor<> *Vijka(
    isArgumentGiven("HHHPCoulombIntegrals") ?
    new Tensor<>(4, ooov, syms, *Cc4s::world, "Vijka") : nullptr
  );
  Tensor<> *Vabci(
    isArgumentGiven("PPPHCoulombIntegrals") ?
    new Tensor<>(4, vvvo, syms, *Cc4s::world, "Vabci") : nullptr
  );
  // Initialization of tensors created from already existing ones
  T *Vijab(
    isArgumentGiven("HHPPCoulombIntegrals") ?
    new Tensor<>(4, oovv, syms, *Cc4s::world, "Vijab") : nullptr
  );
  T *Viajk(
    isArgumentGiven("HPHHCoulombIntegrals") ?
    new Tensor<>(4, ovoo, syms, *Cc4s::world, "Viajk") : nullptr
  );
  T *Vaijb(
    isArgumentGiven("PHHPCoulombIntegrals") ?
    new Tensor<>(4, voov, syms, *Cc4s::world, "Vaijb") : nullptr
  );
  T *Viajb(
    isArgumentGiven("HPHPCoulombIntegrals") ?
    new Tensor<>(4, ovov, syms, *Cc4s::world, "Viajb") : nullptr
  );
  T *Viabc(
    isArgumentGiven("HPPPCoulombIntegrals") ?
    new Tensor<>(4, ovvv, syms, *Cc4s::world, "Viabc") : nullptr
  );
  T *Vabic(
    isArgumentGiven("PPHPCoulombIntegrals") ?
    new Tensor<>(4, vvov, syms, *Cc4s::world, "Vabic") : nullptr
  );

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
  // Allocation of tensors created from already existing ones
  if (Vijab) {
    allocatedTensorArgument("HHPPCoulombIntegrals", Vijab);
  }
  if (Viajk) {
    allocatedTensorArgument("HPHHCoulombIntegrals", Viajk);
  }
  if (Vaijb) {
    allocatedTensorArgument("PHHPCoulombIntegrals", Vaijb);
  }
  if (Viajb) {
    allocatedTensorArgument("HPHPCoulombIntegrals", Viajb);
  }
  if (Viabc) {
    allocatedTensorArgument("HPPPCoulombIntegrals", Viabc);
  }
  if (Vabic) {
    allocatedTensorArgument("PPHPCoulombIntegrals", Vabic);
  }


  // Allocate and compute GammaGab,GammaGai,GammaGij from GammaGqr
  int GaiStart[] = {0 ,aStart,iStart};
  int GaiEnd[]   = {NG,aEnd,  iEnd  };
  int GabStart[] = {0 ,aStart,aStart};
  int GabEnd[]   = {NG,aEnd,  aEnd  };
  int GijStart[] = {0 ,iStart,iStart};
  int GijEnd[]   = {NG,iEnd,  iEnd  };
  Tensor<complex> GammaGai(GammaGqr->slice(GaiStart,GaiEnd));
  Tensor<complex> GammaGab(GammaGqr->slice(GabStart,GabEnd));
  Tensor<complex> GammaGij(GammaGqr->slice(GijStart,GijEnd));

  // Split GammaGab,GammaGai,GammaGia,GammaGij into real and imaginary parts
  Tensor<> realGammaGai(
    3, GammaGai.lens, GammaGai.sym, *GammaGai.wrld, "RealGammaGai"
  );
  Tensor<> imagGammaGai(
    3, GammaGai.lens, GammaGai.sym, *GammaGai.wrld, "ImagGammaGai"
  );
  fromComplexTensor(GammaGai, realGammaGai, imagGammaGai);

  Tensor<> realGammaGab(
    3, GammaGab.lens, GammaGab.sym, *GammaGab.wrld, "RealGammaGab"
  );
  Tensor<> imagGammaGab(
    3, GammaGab.lens, GammaGab.sym, *GammaGab.wrld, "ImagGammaGab"
  );
  fromComplexTensor(GammaGab, realGammaGab, imagGammaGab);

  Tensor<> realGammaGij(
    3, GammaGij.lens, GammaGij.sym, *GammaGij.wrld, "RealGammaGij"
  );
  Tensor<> imagGammaGij(
    3, GammaGij.lens, GammaGij.sym, *GammaGij.wrld, "ImagGammaGij"
  );
  fromComplexTensor(GammaGij, realGammaGij, imagGammaGij);

  // Compute the integrals Vabij Vaibj Vijkl Vabcd
  if (Vabcd) {
    LOG(1, "Integrals") <<
      "Evaluating " << Vabcd->get_name() << std::endl;
    (*Vabcd)["abcd"]  = realGammaGab["Gac"] * realGammaGab["Gbd"];
    (*Vabcd)["abcd"] += imagGammaGab["Gac"] * imagGammaGab["Gbd"];
  }
  if (Vaibj) {
    LOG(1, "Integrals") <<
      "Evaluating " << Vaibj->get_name() << std::endl;
    (*Vaibj)["aibj"]  = realGammaGab["Gab"] * realGammaGij["Gij"];
    (*Vaibj)["aibj"] += imagGammaGab["Gab"] * imagGammaGij["Gij"];
  }
  if (Vabij) {
    LOG(1, "Integrals") <<
      "Evaluating " << Vabij->get_name() << std::endl;
    (*Vabij)["abij"]  = realGammaGai["Gai"] * realGammaGai["Gbj"];
    (*Vabij)["abij"] += imagGammaGai["Gai"] * imagGammaGai["Gbj"];
  }
  if (Vijkl) {
    LOG(1, "Integrals") <<
      "Evaluating " << Vijkl->get_name() << std::endl;
    (*Vijkl)["ijkl"]  = realGammaGij["Gik"] * realGammaGij["Gjl"];
    (*Vijkl)["ijkl"] += imagGammaGij["Gik"] * imagGammaGij["Gjl"];
  }
  if (Vijka) {
    LOG(1, "Integrals") <<
      "Evaluating " << Vijka->get_name() << std::endl;
    (*Vijka)["ijka"]  = realGammaGij["Gik"] * realGammaGai["Gaj"];
    (*Vijka)["ijka"] += imagGammaGij["Gik"] * imagGammaGai["Gaj"];
  }
  if (Vabci) {
    LOG(1, "Integrals") <<
      "Evaluating " << Vabci->get_name() << std::endl;
    (*Vabci)["abci"]  = realGammaGab["Gac"] * realGammaGai["Gbi"];
    (*Vabci)["abci"] += imagGammaGab["Gac"] * imagGammaGai["Gbi"];
  }

  // Create the rest of integrals from the already given ones
  if (Vijab) {
    //                   -----
    // Remember: Vijab = Vabij
    (*Vijab)["ijab"] =  (*Vabij)["abij"];
    if (antisymmetrize) {
      (*Vijab)["ijab"] +=  ( - 1.0) * (*Vabij)["abji"];
    }
    conjugate(*Vijab);
  }
  if (Viajk) {
    //                   -----
    // Remember: Viajk = Vjkia
    (*Viajk)["iajk"] =  (*Vijka)["jkia"];
    if (antisymmetrize) {
      (*Viajk)["iajk"] += ( - 1.0 ) * (*Vijka)["kjia"];
    }
    conjugate(*Viajk);
  }
  // TODO: Do it and check
  if (Viajb) {
    if (antisymmetrize) {
      // FIXME: Assumes real orbitals
      (*Vaijb)["aijb"] = (*Vabij)["abji"];
      (*Viajb)["iajb"] = ( - 1.0 ) * (*Vaijb)["aijb"];
      (*Viajb)["iajb"] +=  (*Vaibj)["aibj"];
    }
  }
  if (Viabc) {
    //                           -----
    // Remember: Viabc = Vaicb = Vcbai
    // TODO: review if it is really like this
    (*Viabc)["iabc"] =  (*Vabci)["abci"];
    conjugate(*Viabc);
    if (antisymmetrize) {
      (*Viabc)["iabc"] -= (*Vabci)["acbi"];
    }
  }
  if (Vabic) {
    // Remember: Vabic = Vbaci
    (*Vabic)["abic"] = (*Vabci)["baci"];
    if (antisymmetrize) {
      (*Vabic)["abic"] = ( - 1.0 ) * (*Vabci)["abci"];
    }
  }

  if (antisymmetrize) {
    // Antisymmetrize integrals calculated directly from Gamma
    // IMPORTANT: This must be written after the creation of the integrals
    //            depending on them.
    if (Vijkl) (*Vijkl)["ijkl"] -= (*Vijkl)["ijlk"];
    if (Vabcd) (*Vabcd)["abcd"] -= (*Vabcd)["abdc"];
    if (Vijka) (*Vijka)["ijka"] -= (*Vijka)["jika"];
    if (Vaibj) (*Vaibj)["aibj"] -= (*Vabij)["baij"];
    if (Vabci) (*Vabci)["abci"] -= (*Vabci)["baci"];
    if (Vabij) (*Vabij)["abij"] -= (*Vabij)["abji"];
  }

  /*
  // debugging info (imaginary part of integrals)
  if (Vabcd) {
    Tensor<> imagVabcd(false, *Vabcd);
    imagVabcd.set_name("imagVabcd");
    LOG(1, "Integrals") << "Evaluating " << imagVabcd.get_name() << std::endl;
    imagVabcd["abcd"]  = realGammaGab["Gac"] * imagGammaGab["Gbd"];
    imagVabcd["abcd"] -= imagGammaGab["Gac"] * realGammaGab["Gbd"];
    double norm(frobeniusNorm(imagVabcd));
    LOG(1, "Integrals") << "Norm of " << imagVabcd.get_name() << " =" << norm << std::endl;
    norm=frobeniusNorm(*Vabcd);
    LOG(1, "Integrals") << "Norm of " << Vabcd->get_name() << " =" << norm << std::endl;
  }

  if (Vaibj) {
    Tensor<> imagVaibj(false, *Vaibj);
    imagVaibj.set_name("imagVaibj");
    LOG(1, "Integrals") << "Evaluating " << imagVaibj.get_name() << std::endl;
    imagVaibj["aibj"]  = realGammaGab["Gab"] * imagGammaGij["Gij"];
    imagVaibj["aibj"] -= imagGammaGab["Gab"] * realGammaGij["Gij"];
    double norm(frobeniusNorm(imagVaibj));
    LOG(1, "Integrals") << "Norm of " << imagVaibj.get_name() << " =" << norm << std::endl;
    norm=frobeniusNorm(*Vaibj);
    LOG(1, "Integrals") << "Norm of " << Vaibj->get_name() << " =" << norm << std::endl;
  }

  if (Vabij) {
    Tensor<> imagVabij(false, *Vabij);
    imagVabij.set_name("imagVabij");
    LOG(1, "Integrals") << "Evaluating " << imagVabij.get_name() << std::endl;
    imagVabij["abij"]  = realGammaGai["Gai"] * imagGammaGai["Gbj"];
    imagVabij["abij"] -= imagGammaGai["Gai"] * realGammaGai["Gbj"];
    double norm(frobeniusNorm(imagVabij));
    LOG(1, "Integrals") << "Norm of " << imagVabij.get_name() << " =" << norm << std::endl;
    norm=frobeniusNorm(*Vabij);
    LOG(1, "Integrals") << "Norm of " << Vabij->get_name() << " =" << norm << std::endl;
  }

  if (Vijkl) {
    Tensor<> imagVijkl(false, *Vijkl);
    imagVijkl.set_name("imagVijkl");
    LOG(1, "Integrals") << "Evaluating " << imagVijkl.get_name() << std::endl;
    imagVijkl["ijkl"]  = realGammaGij["Gik"] * imagGammaGij["Gjl"];
    imagVijkl["ijkl"] -= imagGammaGij["Gik"] * realGammaGij["Gjl"];
    double norm(frobeniusNorm(imagVijkl));
    LOG(1, "Integrals") << "Norm of " << imagVijkl.get_name() << " =" << norm << std::endl;
    norm=frobeniusNorm(*Vijkl);
    LOG(1, "Integrals") << "Norm of " << Vijkl->get_name() << " =" << norm << std::endl;
  }
  */

}

void CoulombIntegralsFromVertex::dryRun() {
  // Read the Coulomb vertex GammaGqr
  DryTensor<complex> *GammaGqr(
    getTensorArgument<complex, DryTensor<complex>>("CoulombVertex")
  );

  // Read the Particle/Hole Eigenenergies
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies")
  );

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

  DryTensor<> *Vabcd(isArgumentGiven("PPPPCoulombIntegrals") ?
    new DryTensor<>(4, vvvv, syms, SOURCE_LOCATION) : nullptr);
  DryTensor<> *Vaibj(isArgumentGiven("PHPHCoulombIntegrals") ?
    new DryTensor<>(4, vovo, syms, SOURCE_LOCATION) : nullptr);
  DryTensor<> *Vabij(isArgumentGiven("PPHHCoulombIntegrals") ?
    new DryTensor<>(4, vvoo, syms, SOURCE_LOCATION) : nullptr);
  DryTensor<> *Vijkl(isArgumentGiven("HHHHCoulombIntegrals") ?
    new DryTensor<>(4, oooo, syms, SOURCE_LOCATION) : nullptr);
  DryTensor<> *Vijka(isArgumentGiven("HHHPCoulombIntegrals") ?
    new DryTensor<>(4, ooov, syms, SOURCE_LOCATION) : nullptr);
  DryTensor<> *Vabci(isArgumentGiven("PPPHCoulombIntegrals") ?
    new DryTensor<>(4, vvvo, syms, SOURCE_LOCATION) : nullptr);

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

  DryTensor<complex> GammaGai(3, GaiLens, syms, SOURCE_LOCATION);
  DryTensor<complex> GammaGab(3, GabLens, syms, SOURCE_LOCATION);
  DryTensor<complex> GammaGij(3, GijLens, syms, SOURCE_LOCATION);

  // Split GammaGab,GammaGai,GammaGij into real and imaginary parts
  DryTensor<> realGammaGai(3, GaiLens, syms, SOURCE_LOCATION);
  DryTensor<> imagGammaGai(3, GaiLens, syms, SOURCE_LOCATION);

  DryTensor<> realGammaGab(3, GabLens, syms, SOURCE_LOCATION);
  DryTensor<> imagGammaGab(3, GabLens, syms, SOURCE_LOCATION);

  DryTensor<> realGammaGij(3, GijLens, syms, SOURCE_LOCATION);
  DryTensor<> imagGammaGij(3, GijLens, syms, SOURCE_LOCATION);
}

