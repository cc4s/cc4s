#include <algorithms/CcsdDiagrammaticDecomposition.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <array>
#include <string>
#include <util/SharedPointer.hpp>


using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CcsdDiagrammaticDecomposition);

CcsdDiagrammaticDecomposition::CcsdDiagrammaticDecomposition(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

CcsdDiagrammaticDecomposition::~CcsdDiagrammaticDecomposition() {
}

//////////////////////////////////////////////////////////////////////
// Hirata iteration routine for the CCSD amplitudes Tabij and Tai from
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
//////////////////////////////////////////////////////////////////////

void CcsdDiagrammaticDecomposition::run() {
  // get singles and doubles part of the amplitudes
  auto Tai(getTensorArgument("CcsdSinglesAmplitudes"));
  Tai->set_name("Tai");
  auto Tabij(getTensorArgument("CcsdDoublesAmplitudes"));
  Tabij->set_name("Tabij");

  int singlesIn(getIntegerArgument("ZeroSinglesIn",-1));
  if (singlesIn > 0) {
    (*Tai)["ai"] = 0.;
  }

  auto Rai( new Tensor<>(*Tai));
  Rai->set_name("Rai");
  (*Rai)["ai"] = 0.;
  auto Rabij( new Tensor<>(*Tabij));
  Rabij->set_name("Rabij");
  (*Rabij)["abij"] = 0.;


  // construct energy denominator
  Tensor<> *epsi(getTensorArgument("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument("ParticleEigenEnergies"));
  auto deltaabij(new Tensor<>(*Rabij));
  deltaabij->set_name("deltaabij");
  (*deltaabij)["abij"]  = (*epsi)["i"];
  (*deltaabij)["abij"] += (*epsi)["j"];
  (*deltaabij)["abij"] -= (*epsa)["a"];
  (*deltaabij)["abij"] -= (*epsa)["b"];
  auto deltaai(new Tensor<>(*Rai));
  deltaai->set_name("deltaai");
  (*deltaai)["ai"]  = (*epsi)["i"];
  (*deltaai)["ai"] -= (*epsa)["a"];


  auto Vabij(getTensorArgument("PPHHCoulombIntegrals"));
  auto Vaibj(getTensorArgument("PHPHCoulombIntegrals"));
  auto Vijkl(getTensorArgument("HHHHCoulombIntegrals"));
  auto Vijka(getTensorArgument("HHHPCoulombIntegrals"));

  // Read the Coulomb vertex GammaGqr
  auto GammaGqr( getTensorArgument<complex>("CoulombVertex"));

  // Compute the No,Nv,NG,Np
  int No(Vabij->lens[2]);
  int Nv(Vabij->lens[0]);
  int NG(GammaGqr->lens[0]);
  int Np(GammaGqr->lens[1]);

  // Allocate and compute GammaGab,GammaGai,GammaGij from GammaGqr
  int GaiStart[] = {0 ,No, 0};
  int GaiEnd[]   = {NG,Np,No};
  int GabStart[] = {0 ,No,No};
  int GabEnd[]   = {NG,Np,Np};
  int GijStart[] = {0 , 0, 0};
  int GijEnd[]   = {NG,No,No};
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

  std::array<int,4> syms({{ NS, NS, NS, NS }});
  std::array<int,4> voov({{ Nv, No, No, Nv }});
  std::array<int,2> vv({{ Nv, Nv }});
  std::array<int,2> vo({{ Nv, No }});
  std::array<int,2> oo({{ No, No }});

  //*************************************************************************
  //****************  T2 amplitude equations  *******************************
  //*************************************************************************


  // Intermediates used both by T1 and T2
  Tensor<> Kac(2, vv.data(), syms.data(), *Vabij->wrld, "Kac");
  Tensor<> Kki(2, oo.data(), syms.data(), *Vabij->wrld, "Kki");
  Tensor<> Xabij(*Tabij);
  Xabij.set_name("Xabij");
  Xabij["abij"] += (*Tai)["ai"] * (*Tai)["bj"];

  {
    // Intermediates used for T2 amplitudes
    Tensor<> Lac(2, vv.data(), syms.data(), *Vabij->wrld, "Lac");
    Tensor<> Lki(2, oo.data(), syms.data(), *Vabij->wrld, "Lki");

    // Build Kac
    Kac["ac"]  = (-2.0) * (*Vabij)["cdkl"] * Xabij["adkl"];
    Kac["ac"] += ( 1.0) * (*Vabij)["dckl"] * Xabij["adkl"];

    // Build Lac
    Lac["ac"]  = Kac["ac"];
    Lac["ac"] += ( 2.0) * realGammaGab["Gca"] * realGammaGai["Gdk"] * (*Tai)["dk"];
    Lac["ac"] += ( 2.0) * imagGammaGab["Gca"] * imagGammaGai["Gdk"] * (*Tai)["dk"];
    Lac["ac"] += (-1.0) * realGammaGai["Gck"] * realGammaGab["Gda"] * (*Tai)["dk"];
    Lac["ac"] += (-1.0) * imagGammaGai["Gck"] * imagGammaGab["Gda"] * (*Tai)["dk"];

    // Build Kki
    Kki["ki"]  = ( 2.0) * (*Vabij)["cdkl"] * Xabij["cdil"];
    Kki["ki"] += (-1.0) * (*Vabij)["dckl"] * Xabij["cdil"];

    // Build Lki
    Lki["ki"]  = ( 1.0) *   Kki   ["ki"];
    Lki["ki"] += ( 2.0) * (*Vijka)["klic"] * (*Tai)["cl"];
    Lki["ki"] += (-1.0) * (*Vijka)["lkic"] * (*Tai)["cl"];

    // Contract Lac with T2 Amplitudes
    (*Rabij)["abij"] = ( 1.0) * Lac["ac"] * (*Tabij)["cbij"];
    (*Rabij)["abij"] += (*Rabij)["baji"];
    evaluateEnergy("Lac", *deltaabij, *deltaai, *Rabij, *Rai);
    // Contract Lki with T2 Amplitudes
    (*Rabij)["abij"] = (-1.0) * Lki["ki"] * (*Tabij)["abkj"];
    (*Rabij)["abij"] += (*Rabij)["baji"];
    evaluateEnergy("Lki", *deltaabij, *deltaai, *Rabij, *Rai);
  }

  {
    // Contract Coulomb integrals with T2 amplitudes
    Tensor<> realDressedGammaGai(realGammaGai);
    Tensor<> imagDressedGammaGai(imagGammaGai);
    realDressedGammaGai.set_name("realDressedGammaGai");
    imagDressedGammaGai.set_name("imagDressedGammaGai");

    realDressedGammaGai["Gai"] += (-1.0) * realGammaGij["Gki"] * (*Tai)["ak"];
    imagDressedGammaGai["Gai"] += (-1.0) * imagGammaGij["Gki"] * (*Tai)["ak"];

    (*Rabij)["abij"] = ( 1.0) * realDressedGammaGai["Gai"] * realGammaGab["Gbc"] * (*Tai)["cj"];
    (*Rabij)["abij"] += ( 1.0) * imagDressedGammaGai["Gai"] * imagGammaGab["Gbc"] * (*Tai)["cj"];

    (*Rabij)["abij"] += (-1.0) * (*Vijka)["jika"] * (*Tai)["bk"];
    (*Rabij)["abij"] += (-1.0) * (*Tai)["bk"] * (*Vabij)["acik"] * (*Tai)["cj"];
    (*Rabij)["abij"] += (*Rabij)["baji"];
    evaluateEnergy("singlesContribtion", *deltaabij, *deltaai, *Rabij, *Rai);
  }

  {
    // Build Xakic
    Tensor<> Xakic(4, voov.data(), syms.data(), *Vabij->wrld, "Xakic");

    Tensor<> realDressedGammaGai(realGammaGai);
    Tensor<> imagDressedGammaGai(imagGammaGai);
    realDressedGammaGai.set_name("realDressedGammaGai");
    imagDressedGammaGai.set_name("imagDressedGammaGai");

    realDressedGammaGai["Gai"] += (-1.0) * realGammaGij["Gil"] * (*Tai)["al"];
    imagDressedGammaGai["Gai"] += (-1.0) * imagGammaGij["Gil"] * (*Tai)["al"];

    realDressedGammaGai["Gai"] += ( 1.0) * realGammaGab["Gad"] * (*Tai)["di"];
    imagDressedGammaGai["Gai"] += ( 1.0) * imagGammaGab["Gad"] * (*Tai)["di"];

    Xakic["akic"]  = ( 1.0) * realDressedGammaGai["Gai"] * realGammaGai["Gck"];
    Xakic["akic"] += ( 1.0) * imagDressedGammaGai["Gai"] * imagGammaGai["Gck"];

    (*Rabij)["abij"] = (2.0) * Xakic["akic"] * (*Tabij)["cbkj"];
    (*Rabij)["abij"] += (*Rabij)["baji"];
    evaluateEnergy("ring", *deltaabij, *deltaai, *Rabij, *Rai);

    (*Rabij)["abij"] = (-1.0) * Xakic["akic"] * (*Tabij)["bckj"];
    (*Rabij)["abij"] += (*Rabij)["baji"];
    evaluateEnergy("strange", *deltaabij, *deltaai, *Rabij, *Rai);


    // Intermediate tensor Yabij=T2-2*T1*T1
    Tensor<> Yabij(*Tabij);
    Yabij.set_name("Yabij");
    Yabij["abij"] += ( 2.0) * (*Tai)["ai"] * (*Tai)["bj"];

    Xakic["akic"] = (-0.5) * (*Vabij)["dclk"] *   Yabij ["dail"];
    Xakic["akic"] += ( 1.0) * (*Vabij)["dclk"] * (*Tabij)["adil"];
    Xakic["akic"] += (-0.5) * (*Vabij)["cdlk"] * (*Tabij)["adil"];

    // Contract and Xakic intermediates with T2 amplitudes Tabij
    (*Rabij)["abij"] = ( 2.0) * Xakic["akic"] * (*Tabij)["cbkj"];
    (*Rabij)["abij"] += (*Rabij)["baji"];
    evaluateEnergy("qu.ring", *deltaabij, *deltaai, *Rabij, *Rai);

    (*Rabij)["abij"] = (-1.0) * Xakic["akic"] * (*Tabij)["bckj"];
    (*Rabij)["abij"] += (*Rabij)["baji"];
    evaluateEnergy("qu.strange", *deltaabij, *deltaai, *Rabij, *Rai);

  }

  {
    // Build Xakci
    Tensor<> Xakci(false, *Vaibj);
    Xakci.set_name("Xakci");

    // Construct dressed Coulomb vertex GammaGab and GammaGij
    Tensor<> realDressedGammaGab(realGammaGab);
    Tensor<> imagDressedGammaGab(imagGammaGab);
    realDressedGammaGab.set_name("realDressedGammaGab");
    imagDressedGammaGab.set_name("imagDressedGammaGab");

    Tensor<> realDressedGammaGij(realGammaGij);
    Tensor<> imagDressedGammaGij(imagGammaGij);
    realDressedGammaGij.set_name("realDressedGammaGij");
    imagDressedGammaGij.set_name("imagDressedGammaGij");

    realDressedGammaGab["Gac"] += (-1.0) * realGammaGai["Gcl"] * (*Tai)["al"];
    imagDressedGammaGab["Gac"] += (-1.0) * imagGammaGai["Gcl"] * (*Tai)["al"];

    realDressedGammaGij["Gki"] += ( 1.0) * realGammaGai["Gdk"] * (*Tai)["di"];
    imagDressedGammaGij["Gki"] += ( 1.0) * imagGammaGai["Gdk"] * (*Tai)["di"];

    // Xakci = Vakci - Vlkci * Tal + Vakcd * Tdi - Vcdlk * Tdail
    Xakci["akci"]  = ( 1.0) * realDressedGammaGab["Gac"] * realDressedGammaGij["Gki"];
    Xakci["akci"] += ( 1.0) * imagDressedGammaGab["Gac"] * imagDressedGammaGij["Gki"];

    (*Rabij)["abij"] = (-1.0) * Xakci["akci"] * (*Tabij)["cbkj"];
    (*Rabij)["abij"] += (*Rabij)["baji"];
    evaluateEnergy("ph1", *deltaabij, *deltaai, *Rabij, *Rai);

    (*Rabij)["abij"] = (-1.0) * Xakci["bkci"] * (*Tabij)["ackj"];
    (*Rabij)["abij"] += (*Rabij)["baji"];
    evaluateEnergy("ph1", *deltaabij, *deltaai, *Rabij, *Rai);

    Xakci["akci"] = (-0.5) * (*Vabij)["cdlk"] * (*Tabij)["dail"];
    (*Rabij)["abij"]  = (-1.0) * Xakci["akci"] * (*Tabij)["cbkj"];
    (*Rabij)["abij"] += (-1.0) * Xakci["bkci"] * (*Tabij)["ackj"];
    (*Rabij)["abij"] += (*Rabij)["baji"];
    evaluateEnergy("qu.ph", *deltaabij, *deltaai, *Rabij, *Rai);
  }

  //////////////////////////////////////////////////////////////////////
  // Now add all terms to Rabij that do not need to be symmetrized with
  // the permutation operator
  //////////////////////////////////////////////////////////////////////

  // Add Vabij to Rabij (MP2 term)
  (*Rabij)["abij"] = (*Vabij)["abij"];
  evaluateEnergy("driver", *deltaabij, *deltaai, *Rabij, *Rai);

  {
    // Build Xklij intermediate

    // dressed hh-ladder...
    Tensor<> Xklij(false, *Vijkl);
    Xklij.set_name("Xklij");

    Xklij["klij"]  = (*Vijkl)["klij"];
    Xklij["klij"] += (*Vijka)["klic"] * (*Tai)["cj"];
    Xklij["klij"] += (*Vijka)["lkjc"] * (*Tai)["ci"];

    // Contract Xklij with T2+T1*T1 Amplitudes via Xabij
    (*Rabij)["abij"] =  Xklij["klij"] * Xabij["abkl"];
    evaluateEnergy("hhladder", *deltaabij, *deltaai, *Rabij, *Rai);

    // Construct last term
    Xklij["klij"]  = (*Vabij)["cdkl"] * Xabij["cdij"];

    (*Rabij)["abij"] = Xklij["klij"] * Xabij["abkl"];
    evaluateEnergy("qu.hh.", *deltaabij, *deltaai, *Rabij, *Rai);
  }

  //********************************************************************************
  //*****************************  T1 contributions  *******************************
  //********************************************************************************
  {
    (*Rabij)["abij"] = 0.;
    // Intermediates used for T1 amplitudes
    Tensor<> Kck(2, vo.data(), syms.data(), *Vabij->wrld, "Kck");

    // Contract Kac and Kki with T1 amplitudes
    (*Rai)["ai"] += ( 1.0) * Kac["ac"] * (*Tai)["ci"];
    (*Rai)["ai"] += (-1.0) * Kki["ki"] * (*Tai)["ak"];

    // Build Kck
    Kck["ck"]  = ( 2.0) * (*Vabij)["cdkl"] * (*Tai)["dl"];
    Kck["ck"] += (-1.0) * (*Vabij)["dckl"] * (*Tai)["dl"];

    // Contract all the rest terms with T1 and T2 amplitudes
    (*Rai)["ai"] += ( 2.0) * Kck["ck"] * (*Tabij)["caki"];
    (*Rai)["ai"] += (-1.0) * Kck["ck"] * (*Tabij)["caik"];
    (*Rai)["ai"] += ( 1.0) * (*Tai)["ak"] * Kck["ck"] * (*Tai)["ci"];
    (*Rai)["ai"] += ( 2.0) * (*Vabij)["acik"] * (*Tai)["ck"];
    (*Rai)["ai"] += (-1.0) * (*Vaibj)["ciak"] * (*Tai)["ck"];

    (*Rai)["ai"] += ( 2.0) * realGammaGab["Gca"] * realGammaGai["Gdk"] * Xabij["cdik"];
    (*Rai)["ai"] += ( 2.0) * imagGammaGab["Gca"] * imagGammaGai["Gdk"] * Xabij["cdik"];
    (*Rai)["ai"] += (-1.0) * realGammaGab["Gda"] * realGammaGai["Gck"] * Xabij["cdik"];
    (*Rai)["ai"] += (-1.0) * imagGammaGab["Gda"] * imagGammaGai["Gck"] * Xabij["cdik"];

    (*Rai)["ai"] += (-2.0) * (*Vijka)["klic"] * Xabij["ackl"];
    (*Rai)["ai"] += ( 1.0) * (*Vijka)["lkic"] * Xabij["ackl"];

    evaluateEnergy("singles", *deltaabij, *deltaai, *Rabij, *Rai);
  }
  //*********************************************************************************
  //***************************  PPL Term *******************************************
  //*********************************************************************************

  // Read the integralsSliceSize. If not provided use No
  int integralsSliceSize(getIntegerArgument("integralsSliceSize"));
  
  int numberSlices(int(ceil(double(Nv)/integralsSliceSize)));

  Tensor<> realDressedGammaGab(realGammaGab);
  Tensor<> imagDressedGammaGab(imagGammaGab);
  realDressedGammaGab.set_name("realDressedGammaGab");
  imagDressedGammaGab.set_name("imagDressedGammaGab");
  // Construct dressed Coulomb vertex GammaGab
  realDressedGammaGab["Gab"] += (-1.0) * realGammaGai["Gbk"] * (*Tai)["ak"];
  imagDressedGammaGab["Gab"] += (-1.0) * imagGammaGai["Gbk"] * (*Tai)["ak"];

  std::vector<PTR(CTF::Tensor<double>)> realSlicedGammaGab;
  std::vector<PTR(CTF::Tensor<double>)> imagSlicedGammaGab;
  for (int v(0); v < numberSlices; v++){
    int xStart = v*integralsSliceSize;
    int xEnd = std::min((v+1)*integralsSliceSize,Nv);

    int sliceStart[] = { 0, xStart, 0};
    int sliceEnd[]   = { NG, xEnd, Nv};
    realSlicedGammaGab.push_back(
      NEW(CTF::Tensor<double>,realDressedGammaGab.slice(sliceStart,sliceEnd))
    );
    imagSlicedGammaGab.push_back(
      NEW(CTF::Tensor<double>,imagDressedGammaGab.slice(sliceStart,sliceEnd))
    );
  }

  //slice loop starts here
  for (int m(0); m < numberSlices; m++){
    for (int n(m); n < numberSlices; n++){
      int lenscd[] = {
        realSlicedGammaGab[n]->lens[1], realSlicedGammaGab[m]->lens[1],
        realSlicedGammaGab[n]->lens[2], realSlicedGammaGab[m]->lens[2]
      };
      int syms[] = {NS, NS, NS, NS};
      auto Vxycd(new Tensor<>(4, lenscd, syms, *realDressedGammaGab.wrld, "Vxycd"));
 
      (*Vxycd)["xycd"]  = (*realSlicedGammaGab[n])["Gxc"] * (*realSlicedGammaGab[m])["Gyd"];
      (*Vxycd)["xycd"] += (*imagSlicedGammaGab[n])["Gxc"] * (*imagSlicedGammaGab[m])["Gyd"];
 
      int lensij[] = { Vxycd->lens[0], Vxycd->lens[1], No, No};
      Tensor<> Rxyij(4, lensij, syms, *Vxycd->wrld, "Rxyij");
 
      // Contract sliced Vxycd with T2 and T1 Amplitudes using Xabij
      Rxyij["xyij"] = (*Vxycd)["xycd"] * Xabij["cdij"];
      int a(n*integralsSliceSize);
      int b(m*integralsSliceSize);
 
      sliceIntoResiduum(Rxyij, a, b, *Rabij);
      // The integrals of this slice are not needed anymore
      delete Vxycd;
 
    }
  }

  Tensor<> Xklij(false, *Vijkl);
  Xklij.set_name("Xklij");
  // we have to remove the singles sinlges term which appears in our abcd contraction
  (*Rai)["ai"] = 0.0; 
  Xklij["klij"]  = (*Vabij)["cdkl"] * Xabij["cdij"];
  (*Rabij)["abij"] += (-1.0) * Xklij["klij"] * (*Tai)["ai"] * (*Tai)["bj"];

  evaluateEnergy("Ppl", *deltaabij, *deltaai, *Rabij, *Rai);


}


void CcsdDiagrammaticDecomposition::evaluateEnergy(
  std::string diagramType,
  CTF::Tensor<> &deltaabij, CTF::Tensor<> &deltaai,
  CTF::Tensor<> &Rabij,  CTF::Tensor<> &Rai
){
  // energy denominator
  CTF::Transform<double, double>(
    std::function<void(double, double &)>(
      [](double denominator, double &res){
        res = res / denominator;
      }
     )
  )(
    deltaabij["abij"], Rabij["abij"]
  );
  CTF::Transform<double, double>(
    std::function<void(double, double &)>(
      [](double edenom, double &res){
        res = res / edenom;
      }
    )
   )(
    deltaai["ai"], Rai["ai"]
  );

  // get the Coulomb integrals to compute the energy
  auto Vijab(getTensorArgument<>("HHPPCoulombIntegrals"));
  
  Scalar<> energy(*Vijab->wrld);
  //direct term
  energy[""]  = 2.0 * Rabij["abij"] * (*Vijab)["ijab"];
  energy[""] += 2.0 * Rai["ai"] * Rai["bj"] * (*Vijab)["ijab"];
  double dire(energy.get_val());
  // exchange term
  energy[""] =  ( -1.0 ) * Rabij["abij"] * (*Vijab)["ijba"];
  energy[""] += ( -1.0 ) * Rai["ai"] * Rai["bj"] * (*Vijab)["ijba"];
  double exce(energy.get_val());
  LOG(1, diagramType) << std::setprecision(10) <<
    "sing= " << 0.25*dire - 0.5*exce << std::endl;
  LOG(1, diagramType) << std::setprecision(10) <<
    "trip= " << 0.75*dire + 1.5*exce << std::endl;
  LOG(0, diagramType) << std::setprecision(10) <<
    "energy= " << dire + exce << std::endl;

}

void CcsdDiagrammaticDecomposition::sliceIntoResiduum(
  Tensor<> &Rxyij, int a, int b, Tensor<> &Rabij
) {
  int Nx(Rxyij.lens[0]);
  int Ny(Rxyij.lens[1]);
  int No(Rxyij.lens[2]);
  int dstStart[] = { a, b, 0, 0 };
  int dstEnd[] = { a+Nx, b+Ny, No, No };
  int srcStart[] = { 0, 0, 0, 0 };
  int srcEnd[] = { Nx, Ny, No, No };
  // R["abij"] += R["xyij"] at current x,y
  Rabij.slice(dstStart,dstEnd,1.0, Rxyij,srcStart,srcEnd,1.0);
  if (a>b) {
    // Add the same slice at (b,a,j,i):
    dstStart[0] = b; dstStart[1] = a;
    dstEnd[0] = b+Ny; dstEnd[1] = a+Nx;
    srcEnd[0] = Ny; srcEnd[1] = Nx;
    // Swap xy and ij simultaneously
    Tensor<> Ryxji(4, srcEnd, Rxyij.sym, *Rxyij.wrld, "Ryxji");
    Ryxji["yxji"] = Rxyij["xyij"];
    // Add Ryxij to Rabij
    Rabij.slice(dstStart,dstEnd,1.0, Ryxji,srcStart,srcEnd,1.0);
  }
}

