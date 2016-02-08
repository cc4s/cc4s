#include <CcsdCoulombIntegrals.hpp>
#include <util/Complex.hpp>
#include <util/ComplexTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>

using namespace cc4s;
using namespace CTF;

ALGORITHM_REGISTRAR_DEFINITION(CcsdCoulombIntegrals);

CcsdCoulombIntegrals::CcsdCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

CcsdCoulombIntegrals::~CcsdCoulombIntegrals() {
}

/**
 * \brief Calculates Coulomb integrals Vabcd,Vabij,Vaibj,Vabci,Vijka,Vijkl from GammaGai,GammaGab,GammaGij Coulomb Vertices
 */
void CcsdCoulombIntegrals::run() {
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

  LOG(0) <<
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
  Tensor<> *Vabcd(new Tensor<>(4, vvvv, syms, *Cc4s::world, "Vabcd"));
  Tensor<> *Vaibj(new Tensor<>(4, vovo, syms, *Cc4s::world, "Vaibj"));
  Tensor<> *Vabij(new Tensor<>(4, vvoo, syms, *Cc4s::world, "Vabij"));
  Tensor<> *Vijkl(new Tensor<>(4, oooo, syms, *Cc4s::world, "Vijkl"));
  Tensor<> *Vijka(new Tensor<>(4, ooov, syms, *Cc4s::world, "Vijka"));
  Tensor<> *Vabci(new Tensor<>(4, vvvo, syms, *Cc4s::world, "Vabci"));
  allocatedTensorArgument("PPPPCoulombIntegrals", Vabcd);
  allocatedTensorArgument("PHPHCoulombIntegrals", Vaibj);
  allocatedTensorArgument("PPHHCoulombIntegrals", Vabij);
  allocatedTensorArgument("HHHHCoulombIntegrals", Vijkl);
  allocatedTensorArgument("HHHPCoulombIntegrals", Vijka);
  allocatedTensorArgument("PPPHCoulombIntegrals", Vabci);

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
  (*Vabcd)["abcd"]  = realGammaGab["Gac"] * realGammaGab["Gbd"];
  (*Vabcd)["abcd"] += imagGammaGab["Gac"] * imagGammaGab["Gbd"];
  (*Vaibj)["aibj"]  = realGammaGab["Gab"] * realGammaGij["Gij"];
  (*Vaibj)["aibj"] += imagGammaGab["Gab"] * imagGammaGij["Gij"];
  (*Vabij)["abij"]  = realGammaGai["Gai"] * realGammaGai["Gbj"];
  (*Vabij)["abij"] += imagGammaGai["Gai"] * imagGammaGai["Gbj"];
  (*Vijkl)["ijkl"]  = realGammaGij["Gik"] * realGammaGij["Gjl"];
  (*Vijkl)["ijkl"] += imagGammaGij["Gik"] * imagGammaGij["Gjl"];
  (*Vijka)["ijka"]  = realGammaGij["Gik"] * realGammaGai["Gaj"];
  (*Vijka)["ijka"] += imagGammaGij["Gik"] * imagGammaGai["Gaj"];
  (*Vabci)["abci"]  = realGammaGab["Gac"] * realGammaGai["Gbi"];
  (*Vabci)["abci"] += imagGammaGab["Gac"] * imagGammaGai["Gbi"];

  // Print okay
  LOG(0) << " OK" << std::endl;

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
