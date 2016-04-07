#include <algorithms/CoulombIntegralsFromVertex.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
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
  Tensor<complex> *GammaGpq(
    getTensorArgument<complex>("CoulombVertex")
  );

  // Read the Particle/Hole Eigenenergies
  Tensor<> *epsi(
    getTensorArgument<>("HoleEigenEnergies")
  );
  Tensor<> *epsa(
    getTensorArgument<>("ParticleEigenEnergies")
  );

  LOG(0, "CoulombIntegrals") <<
    "Reading Coulomb integrals form vertex " << GammaGpq->get_name() << " ...";

  // Compute the no,nv,nG,np
  int nG(GammaGpq->lens[0]);
  int no(epsi->lens[0]);
  int nv(epsa->lens[0]);
  int np = no + nv;

  // Allocate coulomb integrals Vabij Vaibj Vaijb Vijkl Vabcd
  int syms[] = { NS, NS, NS, NS };
  int vvvv[] = { nv, nv, nv, nv };
  int vovo[] = { nv, no, nv, no };
  int vvoo[] = { nv, nv, no, no };
  int oooo[] = { no, no, no, no };
  int ooov[] = { no, no, no, nv };
  int vvvo[] = { nv, nv, nv, no };

  Tensor<> *Vabcd(isArgumentGiven("PPPPCoulombIntegrals") ?
		  new Tensor<>(4, vvvv, syms, *Cc4s::world, "Vabcd") : nullptr
		  );
  Tensor<> *Vaibj(isArgumentGiven("PHPHCoulombIntegrals") ?
		  new Tensor<>(4, vovo, syms, *Cc4s::world, "Vaibj") : nullptr
		  );
  Tensor<> *Vabij(isArgumentGiven("PPHHCoulombIntegrals") ?
		  new Tensor<>(4, vvoo, syms, *Cc4s::world, "Vabij") : nullptr
		  );
  Tensor<> *Vijkl(isArgumentGiven("HHHHCoulombIntegrals") ?
		  new Tensor<>(4, oooo, syms, *Cc4s::world, "Vijkl") : nullptr
		  );
  Tensor<> *Vijka(isArgumentGiven("HHHPCoulombIntegrals") ?
		  new Tensor<>(4, ooov, syms, *Cc4s::world, "Vijka") : nullptr
		  );
  Tensor<> *Vabci(isArgumentGiven("PPPHCoulombIntegrals") ?
		  new Tensor<>(4, vvvo, syms, *Cc4s::world, "Vabci") : nullptr
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

  // Allocate and compute GammaGab,GammaGai,GammaGij from GammaGpq
  int GaiStart[] = {0 ,no, 0};
  int GaiEnd[]   = {nG,np,no};
  int GabStart[] = {0 ,no,no};
  int GabEnd[]   = {nG,np,np};
  int GijStart[] = {0 , 0, 0};
  int GijEnd[]   = {nG,no,no};
  Tensor<complex> GammaGai(GammaGpq->slice(GaiStart,GaiEnd));
  Tensor<complex> GammaGab(GammaGpq->slice(GabStart,GabEnd));
  Tensor<complex> GammaGij(GammaGpq->slice(GijStart,GijEnd));

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

  // Compute the integrals Vabij Vaibj Vaijb Vijkl Vabcd
  if (Vabcd) {
    (*Vabcd)["abcd"]  = realGammaGab["Gac"] * realGammaGab["Gbd"];
    (*Vabcd)["abcd"] += imagGammaGab["Gac"] * imagGammaGab["Gbd"];
  }
  if (Vaibj) {
    (*Vaibj)["aibj"]  = realGammaGab["Gab"] * realGammaGij["Gij"];
    (*Vaibj)["aibj"] += imagGammaGab["Gab"] * imagGammaGij["Gij"];
  }
  if (Vabij) {
    (*Vabij)["abij"]  = realGammaGai["Gai"] * realGammaGai["Gbj"];
    (*Vabij)["abij"] += imagGammaGai["Gai"] * imagGammaGai["Gbj"];
  }
  if (Vijkl) {
    (*Vijkl)["ijkl"]  = realGammaGij["Gik"] * realGammaGij["Gjl"];
    (*Vijkl)["ijkl"] += imagGammaGij["Gik"] * imagGammaGij["Gjl"];
  }
  if (Vijka) {
    (*Vijka)["ijka"]  = realGammaGij["Gik"] * realGammaGai["Gaj"];
    (*Vijka)["ijka"] += imagGammaGij["Gik"] * imagGammaGai["Gaj"];
  }
  if (Vabci) {
    (*Vabci)["abci"]  = realGammaGab["Gac"] * realGammaGai["Gbi"];
    (*Vabci)["abci"] += imagGammaGab["Gac"] * imagGammaGai["Gbi"];
  }

  // Print okay
  //LOG(0, "CoulombIntegrals") << " OK" << std::endl;

  // Print test norm2 of GammaGai, GammaGab, GammaGij
  // GammaGai
  //double error(realGammaGai.norm2());
  //LOG(4) << "|realGammaGai| = " << error << std::endl;
  //error = imagGammaGai.norm2();
  //LOG(4) << "|imagGammaGai| = " << error << std::endl;
  // GammaGab
  //error = realGammaGab.norm2();
  //LOG(4) << "|realGammaGab| = " << error << std::endl;
  //error = imagGammaGab.norm2();
  //LOG(4) << "|imagGammaGab| = " << error << std::endl;
  // GammaGij
  //error = realGammaGij.norm2();
  //LOG(4) << "|realGammaGij| = " << error << std::endl;
  //error = imagGammaGij.norm2();
  //LOG(4) << "|imagGammaGij| = " << error << std::endl;

  // Print test norm2 of Vabij Vaibj Vaijb Vijkl Vabcd
  //error = Vabcd->norm2();
  //LOG(4) << "|Vabcd| = " << error << std::endl;
  //error = Vabij->norm2();
  //LOG(4) << "|Vabij| = " << error << std::endl;
  //error = Vaibj->norm2();
  //LOG(4) << "|Vaibj| = " << error << std::endl;
  //error = Vijkl->norm2();
  //LOG(4) << "|Vijkl| = " << error << std::endl;
  //error = Vijka->norm2();
  //LOG(4) << "|Vijka| = " << error << std::endl;
  //error = Vabci->norm2();
  //LOG(4) << "|Vabci| = " << error << std::endl;
}

void CoulombIntegralsFromVertex::dryRun() {
  DryTensor<complex> *GammaGpq(
    getTensorArgument<complex, DryTensor<complex>>("CoulombVertex")
  );

  // Read the Particle/Hole Eigenenergies
  DryTensor<> *epsi(
    getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies")
  );
  DryTensor<> *epsa(
    getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies")
  );

  // Compute the no,nv,nG,np
  int nG(GammaGpq->lens[0]);
  int no(epsi->lens[0]);
  int nv(epsa->lens[0]);
  int np = no + nv;

  // Allocate coulomb integrals Vabij Vaibj Vaijb Vijkl Vabcd
  int syms[] = { NS, NS, NS, NS };
  int vvvv[] = { nv, nv, nv, nv };
  int vovo[] = { nv, no, nv, no };
  int vvoo[] = { nv, nv, no, no };
  int oooo[] = { no, no, no, no };
  int ooov[] = { no, no, no, nv };
  int vvvo[] = { nv, nv, nv, no };

  DryTensor<> *Vabcd(isArgumentGiven("PPPPCoulombIntegrals") ?
		  new DryTensor<>(4, vvvv, syms) : nullptr
		  );
  DryTensor<> *Vaibj(isArgumentGiven("PHPHCoulombIntegrals") ?
		  new DryTensor<>(4, vovo, syms) : nullptr
		  );
  DryTensor<> *Vabij(isArgumentGiven("PPHHCoulombIntegrals") ?
		  new DryTensor<>(4, vvoo, syms) : nullptr
		  );
  DryTensor<> *Vijkl(isArgumentGiven("HHHHCoulombIntegrals") ?
		  new DryTensor<>(4, oooo, syms) : nullptr
		  );
  DryTensor<> *Vijka(isArgumentGiven("HHHPCoulombIntegrals") ?
		  new DryTensor<>(4, ooov, syms) : nullptr
		  );
  DryTensor<> *Vabci(isArgumentGiven("PPPHCoulombIntegrals") ?
		  new DryTensor<>(4, vvvo, syms) : nullptr
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

  // Allocate and compute GammaGab,GammaGai,GammaGij from GammaGpq
  int GaiLens[]   = {nG,nv,no};
  int GabLens[]   = {nG,nv,nv};
  int GijLens[]   = {nG,no,no};

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
