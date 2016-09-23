#include <algorithms/CoulombIntegralsFromVertex.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <util/DryTensor.hpp>
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
  // Read the Coulomb vertex GammaGpq
  Tensor<complex> *GammaGpq( getTensorArgument<complex>("CoulombVertex"));

  // Read the Particle/Hole Eigenenergies
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));

  LOG(0, "Integrals") <<
    "Reading Coulomb integrals form vertex " << GammaGpq->get_name() 
					     << std::endl;

  // Compute the No,Nv,NG,Np
  int NG(GammaGpq->lens[0]);
  int No(epsi->lens[0]);
  int Nv(epsa->lens[0]);
  int Np = No + Nv;

  // Allocate coulomb integrals Vabij Vaibj Vaijb Vijkl Vabcd
  int syms[] = { NS, NS, NS, NS };
  int vvvv[] = { Nv, Nv, Nv, Nv };
  int vovo[] = { Nv, No, Nv, No };
  int vvoo[] = { Nv, Nv, No, No };
  int oooo[] = { No, No, No, No };
  int ooov[] = { No, No, No, Nv };
  int vvvo[] = { Nv, Nv, Nv, No };

  Tensor<> *Vabcd(isArgumentGiven("PPPPCoulombIntegrals") ?
		  new Tensor<>(4, vvvv, syms, *Cc4s::world, "Vabcd") : nullptr);
  Tensor<> *Vaibj(isArgumentGiven("PHPHCoulombIntegrals") ?
		  new Tensor<>(4, vovo, syms, *Cc4s::world, "Vaibj") : nullptr);
  Tensor<> *Vabij(isArgumentGiven("PPHHCoulombIntegrals") ?
		  new Tensor<>(4, vvoo, syms, *Cc4s::world, "Vabij") : nullptr);
  Tensor<> *Vijkl(isArgumentGiven("HHHHCoulombIntegrals") ?
		  new Tensor<>(4, oooo, syms, *Cc4s::world, "Vijkl") : nullptr);
  Tensor<> *Vijka(isArgumentGiven("HHHPCoulombIntegrals") ?
		  new Tensor<>(4, ooov, syms, *Cc4s::world, "Vijka") : nullptr);
  Tensor<> *Vabci(isArgumentGiven("PPPHCoulombIntegrals") ?
		  new Tensor<>(4, vvvo, syms, *Cc4s::world, "Vabci") : nullptr);

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

  // Allocate and compute GammaGab,GammaGai,GammaGij from GammaGpq
  int GaiStart[] = {0 ,No, 0};
  int GaiEnd[]   = {NG,Np,No};
  int GabStart[] = {0 ,No,No};
  int GabEnd[]   = {NG,Np,Np};
  int GijStart[] = {0 , 0, 0};
  int GijEnd[]   = {NG,No,No};
  Tensor<complex> GammaGai(GammaGpq->slice(GaiStart,GaiEnd));
  Tensor<complex> GammaGab(GammaGpq->slice(GabStart,GabEnd));
  Tensor<complex> GammaGij(GammaGpq->slice(GijStart,GijEnd));

  // Split GammaGab,GammaGai,GammaGia,GammaGij into real and imaginary parts
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

  Tensor<> realGammaGij(3, GammaGij.lens, GammaGij.sym, 
			*GammaGij.wrld, "RealGammaGij");
  Tensor<> imagGammaGij(3, GammaGij.lens, GammaGij.sym, 
			*GammaGij.wrld, "ImagGammaGij");
  fromComplexTensor(GammaGij, realGammaGij, imagGammaGij);

  // Compute the integrals Vabij Vaibj Vaijb Vijkl Vabcd
  if (Vabcd) {
    LOG(1, "Integrals") << "Evaluating " 
			       << Vabcd->get_name() << std::endl;
    (*Vabcd)["abcd"]  = realGammaGab["Gac"] * realGammaGab["Gbd"];
    (*Vabcd)["abcd"] += imagGammaGab["Gac"] * imagGammaGab["Gbd"];
  }
  if (Vaibj) {
    LOG(1, "Integrals") << "Evaluating " 
			       << Vaibj->get_name() << std::endl;
    (*Vaibj)["aibj"]  = realGammaGab["Gab"] * realGammaGij["Gij"];
    (*Vaibj)["aibj"] += imagGammaGab["Gab"] * imagGammaGij["Gij"];
  }
  if (Vabij) {
    LOG(1, "Integrals") << "Evaluating " 
			       << Vabij->get_name() << std::endl;
    (*Vabij)["abij"]  = realGammaGai["Gai"] * realGammaGai["Gbj"];
    (*Vabij)["abij"] += imagGammaGai["Gai"] * imagGammaGai["Gbj"];
  }
  if (Vijkl) {
    LOG(1, "Integrals") << "Evaluating " 
			       << Vijkl->get_name() << std::endl;
    (*Vijkl)["ijkl"]  = realGammaGij["Gik"] * realGammaGij["Gjl"];
    (*Vijkl)["ijkl"] += imagGammaGij["Gik"] * imagGammaGij["Gjl"];
  }
  if (Vijka) {
    LOG(1, "Integrals") << "Evaluating " 
			       << Vijka->get_name() << std::endl;
    (*Vijka)["ijka"]  = realGammaGij["Gik"] * realGammaGai["Gaj"];
    (*Vijka)["ijka"] += imagGammaGij["Gik"] * imagGammaGai["Gaj"];
  }
  if (Vabci) {
    LOG(1, "Integrals") << "Evaluating " 
			       << Vabci->get_name() << std::endl;
    (*Vabci)["abci"]  = realGammaGab["Gac"] * realGammaGai["Gbi"];
    (*Vabci)["abci"] += imagGammaGab["Gac"] * imagGammaGai["Gbi"];
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
  // Read the Coulomb vertex GammaGpq
  DryTensor<complex> *GammaGpq(getTensorArgument<complex, 
			       DryTensor<complex>>("CoulombVertex"));

  // Read the Particle/Hole Eigenenergies
  DryTensor<> *epsi(getTensorArgument
		    <double, DryTensor<double>>("HoleEigenEnergies"));
  DryTensor<> *epsa(getTensorArgument
		    <double, DryTensor<double>>("ParticleEigenEnergies"));

  // Compute the No,Nv,NG
  int NG(GammaGpq->lens[0]);
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

  // Allocate and compute GammaGab,GammaGai,GammaGij from GammaGpq
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
