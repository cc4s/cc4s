#include <algorithms/CcsdEnergyFromCoulombIntegrals.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <array>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CcsdEnergyFromCoulombIntegrals);

CcsdEnergyFromCoulombIntegrals::CcsdEnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): ClusterSinglesDoublesAlgorithm(argumentList) {
}

CcsdEnergyFromCoulombIntegrals::~CcsdEnergyFromCoulombIntegrals() {
}

void CcsdEnergyFromCoulombIntegrals::iterate(
  int i, Mixer<double> *TaiMixer, Mixer<double> *TabijMixer
) {
  iterate<double>(i, TaiMixer, TabijMixer);
}

void CcsdEnergyFromCoulombIntegrals::iterate(
  int i, Mixer<complex> *TaiMixer, Mixer<complex> *TabijMixer
) {
  iterate<complex>(i, TaiMixer, TabijMixer);
}

//////////////////////////////////////////////////////////////////////
// Hirata iteration routine for the CCSD amplitudes Tabij and Tai from
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
//////////////////////////////////////////////////////////////////////
template <typename F>
void CcsdEnergyFromCoulombIntegrals::iterate(
  int i, Mixer<F> *TaiMixer, Mixer<F> *TabijMixer
) {
  Tensor<F> *Tabij(&TabijMixer->getNext());
  Tabij->set_name("Tabij");
  Tensor<F> *Tai(&TaiMixer->getNext());
  Tai->set_name("Tai");
  // Read all required integrals
  Tensor<F> *Vabij(getTensorArgument<F>("PPHHCoulombIntegrals"));
  Tensor<F> *Vaijb(getTensorArgument<F>("PHHPCoulombIntegrals"));
  Tensor<F> *Vijab(getTensorArgument<F>("HHPPCoulombIntegrals"));
  Tensor<F> *Vaibj(getTensorArgument<F>("PHPHCoulombIntegrals"));
  Tensor<F> *Vijkl(getTensorArgument<F>("HHHHCoulombIntegrals"));
  Tensor<F> *Vijka(getTensorArgument<F>("HHHPCoulombIntegrals"));
  Tensor<F> *Vaijk(getTensorArgument<F>("PHHHCoulombIntegrals"));

  // Read the Coulomb vertex GammaGqr
  Tensor<complex> *GammaGqr( getTensorArgument<complex>("CoulombVertex"));

  // Allocate Tensors for T2 amplitudes
  Tensor<F> Rabij(false, *Tabij);
  Rabij.set_name("Rabij");
  // Allocate Tensors for T1 amplitudes
  Tensor<F> Rai(false, *Tai);
  Rai.set_name("Rai");

  std::string abbreviation(getAbbreviation());
  std::transform(abbreviation.begin(), abbreviation.end(), 
                 abbreviation.begin(), ::toupper);

  // Compute the No,Nv,NG,Np
  int NG(GammaGqr->lens[0]);
  int No(Vabij->lens[2]);
  int Nv(Vabij->lens[0]);
  int Np(GammaGqr->lens[1]);

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
  auto GammaGij( new Tensor<complex>(GammaGqr->slice(GijStart, GijEnd)) );
  auto GammaGia( new Tensor<complex>(GammaGqr->slice(GiaStart, GiaEnd)) );
  auto GammaGai( new Tensor<complex>(GammaGqr->slice(GaiStart, GaiEnd)) );
  auto GammaGab( new Tensor<complex>(GammaGqr->slice(GabStart, GabEnd)) );

  Univar_Function<complex> fConj(conj<complex>);

  Tensor<complex> conjTransposeGammaGai(false, *GammaGai);
  conjTransposeGammaGai.sum(1.0,*GammaGia,"Gia", 0.0,"Gai", fConj);
  Tensor<complex> conjTransposeGammaGia(false, *GammaGia);
  conjTransposeGammaGia.sum(1.0,*GammaGai,"Gai", 0.0,"Gia", fConj);
  Tensor<complex> conjTransposeGammaGab(false, *GammaGab);
  conjTransposeGammaGab.sum(1.0,*GammaGab,"Gba", 0.0,"Gab", fConj);
  Tensor<complex> conjTransposeGammaGij(false, *GammaGij);
  conjTransposeGammaGij.sum(1.0,*GammaGij,"Gji", 0.0,"Gij", fConj);

/*
  syms( std::array<int,4>{{ NS, NS, NS, NS }};
  vvvv = std::array<int,4>{{ Nv, Nv, Nv, Nv }};
  vovo = std::array<int,4>{{ Nv, No, Nv, No }};
  vvoo = std::array<int,4>{{ Nv, Nv, No, No }};
  oovv = std::array<int,4>{{ No, No, Nv, Nv }};
  oooo = std::array<int,4>{{ No, No, No, No }};
  ooov = std::array<int,4>{{ No, No, No, Nv }};
  vvvo = std::array<int,4>{{ Nv, Nv, Nv, No }};
*/
  std::array<int,4> syms({{ NS, NS, NS, NS }});
  std::array<int,4> voov({{ Nv, No, No, Nv }});
  std::array<int,2> vv({{ Nv, Nv }});
  std::array<int,2> vo({{ Nv, No }});
  std::array<int,2> oo({{ No, No }});

  // Intermediates used both by T1 and T2
  Tensor<F> Kac(2, vv.data(), syms.data(), *Vabij->wrld, "Kac");
  Tensor<F> Kki(2, oo.data(), syms.data(), *Vabij->wrld, "Kki");

  //********************************************************************************
  //***********************  T2 amplitude equations  *******************************
  //********************************************************************************

  LOG(1, abbreviation) << "Solving T2 Amplitude Equations" << std::endl;

  if (i == 0 && !isArgumentGiven("startingDoublesAmplitudes") ) {
    // For first iteration compute only the MP2 amplitudes 
    // Since Tabij = 0, Vabij is the only non-zero term
    Rabij["abij"] = (*Vabij)["abij"];
      } 
  else {
    // For the rest iterations compute the CCSD amplitudes

    // Intermediate tensor Xabij=T2+T1*T1
    Tensor<F> Xabij(Tabij);
    Xabij.set_name("Xabij");
    Xabij["abij"] += (*Tai)["ai"] * (*Tai)["bj"];

    {
      // Intermediates used for T2 amplitudes
      Tensor<F> Lac(2, vv.data(), syms.data(), *Vabij->wrld, "Lac");
      Tensor<F> Lki(2, oo.data(), syms.data(), *Vabij->wrld, "Lki");

      Tensor<F> Xklij(false, *Vijkl);
      Xklij.set_name("Xklij");
      Tensor<F> Xakci(false, *Vaibj);
      Xakci.set_name("Xakci");

      // Build Kac
      Kac["ac"]  = (-2.0) * (*Vijab)["klcd"] * Xabij["adkl"];
      Kac["ac"] += ( 1.0) * (*Vijab)["kldc"] * Xabij["adkl"];

      // Build Lac
      Lac["ac"]  = Kac["ac"];
      Lac["ac"] += ( 2.0) * conjTransposeGammaGab["Gca"] * (*GammaGai)["Gdk"] * (*Tai)["dk"];
      Lac["ac"] += (-1.0) * conjTransposeGammaGai["Gck"] * (*GammaGab)["Gda"] * (*Tai)["dk"];

      // Build Kki
      Kki["ki"]  = ( 2.0) * (*Vijab)["klcd"] * Xabij["cdil"];
      Kki["ki"] += (-1.0) * (*Vijab)["kldc"] * Xabij["cdil"];

      // Build Lki
      Lki["ki"]  = ( 1.0) *   Kki   ["ki"];
      Lki["ki"] += ( 2.0) * (*Vijka)["klic"] * (*Tai)["cl"];
      Lki["ki"] += (-1.0) * (*Vijka)["lkic"] * (*Tai)["cl"];
    
      // Contract Lac with T2 Amplitudes
      Rabij["abij"]  = ( 1.0) * Lac["ac"] * (*Tabij)["cbij"];

      // Contract Lki with T2 Amplitudes
      Rabij["abij"] += (-1.0) * Lki["ki"] * (*Tabij)["abkj"];

      {
        // Contract Coulomb integrals with T2 amplitudes
        Tensor<complex> conjTransposeDressedGammaGai(conjTransposeGammaGai);
        conjTransposeDressedGammaGai.set_name("conjTransposeDressedGammaGai");
        conjTransposeDressedGammaGai["Gai"] += (-1.0) * conjTransposeGammaGij["Gki"] * (*Tai)["ak"];
        Rabij["abij"] += ( 1.0) * conjTransposeDressedGammaGai["Gai"] * (*GammaGab)["Gbc"] * (*Tai)["cj"];

        Rabij["abij"] += (-1.0) * (*Vaijk)["akij"] * (*Tai)["bk"];
        Rabij["abij"] += (-1.0) * (*Tai)["bk"] * (*Vaijb)["akic"] * (*Tai)["cj"];
      }
        
      {
        // Build Xakic
        Tensor<F> Xakic(4, voov.data(), syms.data(), *Vabij->wrld, "Xakic");

        Tensor<complex> conjTransposeDressedGammaGai(conjTransposeGammaGai);
        conjTransposeDressedGammaGai.set_name("conjTransposeDressedGammaGai");

        conjTransposeDressedGammaGai["Gai"] += (-1.0) * conjTransposeGammaGij["Gli"] * (*Tai)["al"];

        conjTransposeDressedGammaGai["Gai"] += ( 1.0) * conjTransposeGammaGab["Gad"] * (*Tai)["di"];

        // Intermediate tensor Yabij=T2-2*T1*T1
        Tensor<F> Yabij(Tabij);
        Yabij.set_name("Yabij");
        Yabij["abij"] += ( 2.0) * (*Tai)["ai"] * (*Tai)["bj"];

        conjTransposeDressedGammaGai["Gai"] += (-0.5) * conjTransposeGammaGia["Gld"] * Yabij["dail"];

        Xakic["akic"]  = ( 1.0) * conjTransposeDressedGammaGai["Gai"] * (*GammaGia)["Gkc"];
        

        Yabij["dclk"]  = ( 1.0) * (*Vijab)["lkdc"];
        Yabij["dclk"] += (-0.5) * (*Vijab)["lkcd"];
        
        Xakic["akic"] += Yabij["dclk"] * (*Tabij)["adil"];

        Yabij["cbkj"]  = ( 2.0) * (*Tabij)["cbkj"];
        Yabij["cbkj"] += (-1.0) * (*Tabij)["bckj"];

        Rabij["abij"] += ( 1.0) * Xakic["akic"] * Yabij["cbkj"];
      }
      

      {
        // Build Xakci
        // Construct dressed Coulomb vertex GammaGab and GammaGij
        Tensor<complex> conjTransposeDressedGammaGab(conjTransposeGammaGab);
        conjTransposeDressedGammaGab.set_name("conjTransposeDressedGammaGab");

        Tensor<complex> DressedGammaGij(GammaGij);
        DressedGammaGij.set_name("DressedGammaGij");
        
        conjTransposeDressedGammaGab["Gac"] += (-1.0) * conjTransposeGammaGia["Glc"] * (*Tai)["al"];

        DressedGammaGij["Gki"] += ( 1.0) * (*GammaGia)["Gkd"] * (*Tai)["di"];
            
        // Xakci = Vakci - Vlkci * Tal + Vakcd * Tdi - Vlkcd * Tdail
        Xakci["akci"]  = ( 1.0) * conjTransposeDressedGammaGab["Gac"] * DressedGammaGij["Gki"];

        // Xakci = 0.5 * Vlkcd * Tdail
        Xakci["akci"] += (-0.5) * (*Vijab)["lkcd"] * (*Tabij)["dail"];

        Rabij["abij"] += (-1.0) * Xakci["akci"] * (*Tabij)["cbkj"];
        Rabij["abij"] += (-1.0) * Xakci["bkci"] * (*Tabij)["ackj"];

        // Symmetrize Rabij by applying permutation operator
        Xakci["aibj"]  = Rabij["abij"];
        Rabij["abij"] += Xakci["bjai"]; 
      }

      //////////////////////////////////////////////////////////////////////
      // Now add all terms to Rabij that do not need to be symmetrized with
      // the permutation operator
      //////////////////////////////////////////////////////////////////////

      // Add Vabij to Rabij (MP2 term)
      Rabij["abij"] += (*Vabij)["abij"];

      {
        // Build Xklij intermediate
        Xklij["klij"]  = (*Vijkl)["klij"];
        Xklij["klij"] += (*Vijka)["klic"] * (*Tai)["cj"];
        Xklij["klij"] += (*Vijka)["lkjc"] * (*Tai)["ci"];

        // Contract Xklij with T2+T1*T1 Amplitudes via Xabij
        Rabij["abij"] +=  Xklij["klij"] * Xabij["abkl"];

        // Construct last term
        Xklij["klij"]  = (*Vijab)["klcd"] * Xabij["cdij"];

        // Add last term contracted only with the doubles
        // The singles term is computed in the slicing
        Rabij["abij"] +=  Xklij["klij"] * (*Tabij)["abkl"];
      }


      if (isArgumentGiven("CoulombFactors")) {

	// Read the factorsSliceSize.
	Tensor<complex> *LambdaGR(getTensorArgument<complex>("CoulombFactors"));
	LambdaGR->set_name("LambdaGR");

	int NR(LambdaGR->lens[1]);

	//	factorsSliceSize = getIntegerArgument("factorsSliceSize", NR);

	// calculate decomposition rank
	// if rank is not given use rank factors (if they are not given use rankFactors=2.0)
	//	if (factorsSliceSize == -1) {
	//	  double factorsSliceFactor(getRealArgument("factorsSliceFactor", 1.0));
	//	  rank = NG * rankFactor;
	//	}
                
	int factorsSliceSize(
	  getIntegerArgument("factorsSliceSize", NR)
	);

	// Slice loop starts here
	for (int b(0); b < NR; b += factorsSliceSize) {
	  for (int a(0); a < NR; a += factorsSliceSize) {
	    LOG(1, abbreviation) << "Evaluting Fabij at R=" << a << ", S=" << b << std::endl;
	    Tensor<F> *Fabij(
              sliceAmplitudesFromCoupledCoulombFactors(
              TaiMixer, TabijMixer, a, b, factorsSliceSize
              )
            );
	    Fabij->set_name("Fabij");
	    Rabij["abij"] += (*Fabij)["abij"];
	    delete Fabij;
	  }
	}
      }
      else {
        // Read the integralsSliceSize. If not provided use No
        int integralsSliceSize(getIntegerArgument("integralsSliceSize",No));

        // Slice loop starts here
        for (int b(0); b < Nv; b += integralsSliceSize) {
          for (int a(b); a < Nv; a += integralsSliceSize) {
            LOG(1, abbreviation) << "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
            Tensor<F> *Vxycd(
              sliceCoupledCoulombIntegrals(TaiMixer, a, b, integralsSliceSize)
            );
            Vxycd->set_name("Vxycd");
            int lens[] = { Vxycd->lens[0], Vxycd->lens[1], No, No };
            int syms[] = {NS, NS, NS, NS};
            Tensor<F> Rxyij(4, lens, syms, *Vxycd->wrld, "Rxyij");

            // Contract sliced Vxycd with T2 and T1 Amplitudes using Xabij
            Rxyij["xyij"] = (*Vxycd)["xycd"] * Xabij["cdij"];

            sliceIntoResiduum(Rxyij, a, b, Rabij);
            // The integrals of this slice are not needed anymore
            delete Vxycd;
          }
        }
      }
    }
  }
  
  //********************************************************************************
  //***********************  T1 amplitude equations  *******************************
  //********************************************************************************
  {
    LOG(1, abbreviation) << "Solving T1 Amplitude Equations" << std::endl;
    {  
      // Intermediates used for T1 amplitudes
      Tensor<F> Kck(2, vo.data(), syms.data(), *Vabij->wrld, "Kck");

      // Intermediate tensor Xabij=T2+T1*T1
      Tensor<F> Xabij(Tabij);
      Xabij.set_name("Xabij");
      Xabij["abij"] += (*Tai)["ai"] * (*Tai)["bj"];

      // Contract Kac and Kki with T1 amplitudes
      Rai["ai"]  = ( 1.0) * Kac["ac"] * (*Tai)["ci"];
      Rai["ai"] += (-1.0) * Kki["ki"] * (*Tai)["ak"];

      // Build Kck
      Kck["ck"]  = ( 2.0) * (*Vijab)["klcd"] * (*Tai)["dl"];
      Kck["ck"] += (-1.0) * (*Vijab)["kldc"] * (*Tai)["dl"];

      // Contract all the rest terms with T1 and T2 amplitudes
      Rai["ai"] += ( 2.0) * Kck["ck"] * (*Tabij)["caki"];
      Rai["ai"] += (-1.0) * Kck["ck"] * (*Tabij)["caik"];
      Rai["ai"] += ( 1.0) * (*Tai)["ak"] * Kck["ck"] * (*Tai)["ci"];
      Rai["ai"] += ( 2.0) * (*Vaijb)["akic"] * (*Tai)["ck"];
      Rai["ai"] += (-1.0) * (*Vaibj)["akci"] * (*Tai)["ck"];

      Rai["ai"] += ( 2.0) * conjTransposeGammaGab["Gac"] * (*GammaGia)["Gkd"] * Xabij["cdik"];
      Rai["ai"] += (-1.0) * conjTransposeGammaGab["Gad"] * (*GammaGia)["Gkc"] * Xabij["cdik"];

      Rai["ai"] += (-2.0) * (*Vijka)["klic"] * Xabij["ackl"];
      Rai["ai"] += ( 1.0) * (*Vijka)["lkic"] * Xabij["ackl"];
    }
  }
  
  // Calculate the amplitudes from the residuum
  amplitudesFromResiduum(Rabij, "abij");
  amplitudesFromResiduum(Rai, "ai");
  // Append amplitudes to the mixer
  TabijMixer->append(Rabij);
  TaiMixer->append(Rai);
}

