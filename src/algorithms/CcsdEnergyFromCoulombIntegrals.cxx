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

//////////////////////////////////////////////////////////////////////
// Hirata iteration routine for the CCSD amplitudes Tabij and Tai from
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
//////////////////////////////////////////////////////////////////////

PTR(FockVector<double>) CcsdEnergyFromCoulombIntegrals::getResiduum(
  const int i, const PTR(const FockVector<double>) &amplitudes
) {
  // get singles and doubles part of the amplitudes
  auto Tai( amplitudes->get(0) );
  Tai->set_name("Tai");
  auto Tabij( amplitudes->get(1) );
  Tabij->set_name("Tabij");

  // create residuum and get their singles and doubles part
  auto residuum( NEW(FockVector<double>, *amplitudes) );
  *residuum *= 0.0;
  auto Rai( residuum->get(0) );
  Rai->set_name("Rai");
  auto Rabij( residuum->get(1) );
  Rabij->set_name("Rabij");

  // get part of Coulomb integrals used whether the amplitudes are zero or not
  auto Vabij(getTensorArgument("PPHHCoulombIntegrals"));

  if (i == 0 && !isArgumentGiven("initialDoublesAmplitudes"))  {
    // For first iteration compute only the MP2 amplitudes
    // Since Tabij = 0, Vabij is the only non-zero term
    LOG(1, getCapitalizedAbbreviation()) << "MP2 T2 Amplitudes" << std::endl;
    (*Rabij)["abij"] = (*Vabij)["abij"];
  } else {
    // For the rest iterations compute the CCSD amplitudes

    // Read all required integrals
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

    int distinguishable(
      getIntegerArgument("distinguishable", DEFAULT_DISTINGUISHABLE)
    );

    //*************************************************************************
    //****************  T2 amplitude equations  *******************************
    //*************************************************************************

    LOG(1, getCapitalizedAbbreviation()) <<
      "Solving T2 Amplitude Equations" << std::endl;

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
      if (!distinguishable) {
        Lac["ac"]  = Kac["ac"];
      } else {
        // Intermediate tensor Yabij=T2-2*T1*T1
        Tensor<> Yabij(*Tabij);
        Yabij.set_name("Yabij");
        Yabij["abij"] += ( 2.0) * (*Tai)["ai"] * (*Tai)["bj"];
        Lac["ac"]      = (-1.0) * (*Vabij)["cdkl"] * Yabij["adkl"]; // Use Yabij in DCSD
        Lac["ac"]     += ( 0.5) * (*Vabij)["dckl"] * Yabij["adkl"]; // Use Yabij in DCSD
      }
      Lac["ac"] += ( 2.0) * realGammaGab["Gca"] * realGammaGai["Gdk"] * (*Tai)["dk"];
      Lac["ac"] += ( 2.0) * imagGammaGab["Gca"] * imagGammaGai["Gdk"] * (*Tai)["dk"];
      Lac["ac"] += (-1.0) * realGammaGai["Gck"] * realGammaGab["Gda"] * (*Tai)["dk"];
      Lac["ac"] += (-1.0) * imagGammaGai["Gck"] * imagGammaGab["Gda"] * (*Tai)["dk"];

      // Build Kki
      Kki["ki"]  = ( 2.0) * (*Vabij)["cdkl"] * Xabij["cdil"];
      Kki["ki"] += (-1.0) * (*Vabij)["dckl"] * Xabij["cdil"];

      // Build Lki
      if (!distinguishable) {
        Lki["ki"]  = ( 1.0) *   Kki   ["ki"];
      } else {
        // Intermediate tensor Yabij=T2-2*T1*T1
        Tensor<> Yabij(*Tabij);
        Yabij.set_name("Yabij");
        Yabij["abij"] += ( 2.0) * (*Tai)["ai"] * (*Tai)["bj"];
        Lki["ki"]      = ( 1.0) * (*Vabij)["cdkl"] * Yabij["cdil"]; // Use Yabij in DCSD
        Lki["ki"]     += (-0.5) * (*Vabij)["dckl"] * Yabij["cdil"]; // Use Yabij in DCSD
      }
      Lki["ki"] += ( 2.0) * (*Vijka)["klic"] * (*Tai)["cl"];
      Lki["ki"] += (-1.0) * (*Vijka)["lkic"] * (*Tai)["cl"];

      // Contract Lac with T2 Amplitudes
      (*Rabij)["abij"] += ( 1.0) * Lac["ac"] * (*Tabij)["cbij"];

      // Contract Lki with T2 Amplitudes
      (*Rabij)["abij"] += (-1.0) * Lki["ki"] * (*Tabij)["abkj"];
    }

    {
      // Contract Coulomb integrals with T2 amplitudes
      Tensor<> realDressedGammaGai(realGammaGai);
      Tensor<> imagDressedGammaGai(imagGammaGai);
      realDressedGammaGai.set_name("realDressedGammaGai");
      imagDressedGammaGai.set_name("imagDressedGammaGai");


      realDressedGammaGai["Gai"] += (-1.0) * realGammaGij["Gki"] * (*Tai)["ak"];
      imagDressedGammaGai["Gai"] += (-1.0) * imagGammaGij["Gki"] * (*Tai)["ak"];

      (*Rabij)["abij"] += ( 1.0) * realDressedGammaGai["Gai"] * realGammaGab["Gbc"] * (*Tai)["cj"];
      (*Rabij)["abij"] += ( 1.0) * imagDressedGammaGai["Gai"] * imagGammaGab["Gbc"] * (*Tai)["cj"];

      (*Rabij)["abij"] += (-1.0) * (*Vijka)["jika"] * (*Tai)["bk"];
      (*Rabij)["abij"] += (-1.0) * (*Tai)["bk"] * (*Vabij)["acik"] * (*Tai)["cj"];
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

      // FIXME: there is a better way for the contractions (see complex code)
      Xakic["akic"]  = ( 1.0) * realDressedGammaGai["Gai"] * realGammaGai["Gck"];
      Xakic["akic"] += ( 1.0) * imagDressedGammaGai["Gai"] * imagGammaGai["Gck"];

      // Intermediate tensor Yabij=T2-2*T1*T1
      Tensor<> Yabij(*Tabij);
      Yabij.set_name("Yabij");
      Yabij["abij"] += ( 2.0) * (*Tai)["ai"] * (*Tai)["bj"];

      Xakic["akic"] += (-0.5) * (*Vabij)["dclk"] *   Yabij ["dail"];
      Xakic["akic"] += ( 1.0) * (*Vabij)["dclk"] * (*Tabij)["adil"];
      if (!distinguishable) {
        Xakic["akic"] += (-0.5) * (*Vabij)["cdlk"] * (*Tabij)["adil"];
      }
      // Contract and Xakic intermediates with T2 amplitudes Tabij
      Yabij["cbkj"]  = ( 2.0) * (*Tabij)["cbkj"];
      Yabij["cbkj"] += (-1.0) * (*Tabij)["bckj"];

      (*Rabij)["abij"] += ( 1.0) * Xakic["akic"] * Yabij["cbkj"];
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

      // Xakci = 0.5 * Vcdlk * Tdail
      if (!distinguishable) {
        Xakci["akci"] += (-0.5) * (*Vabij)["cdlk"] * (*Tabij)["dail"];
      }

      // Contract and Xakci intermediates with T2 amplitudes Tabij
      (*Rabij)["abij"] += (-1.0) * Xakci["akci"] * (*Tabij)["cbkj"];
      (*Rabij)["abij"] += (-1.0) * Xakci["bkci"] * (*Tabij)["ackj"];

      // Symmetrize Rabij by applying permutation operator
      (*Rabij)["abij"] += (*Rabij)["baji"];
    }

    //////////////////////////////////////////////////////////////////////
    // Now add all terms to Rabij that do not need to be symmetrized with
    // the permutation operator
    //////////////////////////////////////////////////////////////////////

    // Add Vabij to Rabij (MP2 term)
    (*Rabij)["abij"] += (*Vabij)["abij"];

    {
      // Build Xklij intermediate
      Tensor<> Xklij(false, *Vijkl);
      Xklij.set_name("Xklij");

      Xklij["klij"]  = (*Vijkl)["klij"];
      Xklij["klij"] += (*Vijka)["klic"] * (*Tai)["cj"];
      Xklij["klij"] += (*Vijka)["lkjc"] * (*Tai)["ci"];

      // Contract Xklij with T2+T1*T1 Amplitudes via Xabij
      (*Rabij)["abij"] +=  Xklij["klij"] * Xabij["abkl"];

      // Construct last term
      if (!distinguishable) {
        Xklij["klij"]  = (*Vabij)["cdkl"] * Xabij["cdij"];
      } else {
        Xklij["klij"]  = (*Tai)["ci"] * (*Vabij)["cdkl"] * (*Tai)["dj"];
      }

      // Add last term contracted only with the doubles
      // The singles term is computed in the slicing
      if(isArgumentGiven("PPL")){
        (*Rabij)["abij"] += Xklij["klij"] * Xabij["abkl"];
      }
      else{
        (*Rabij)["abij"] +=  Xklij["klij"] * (*Tabij)["abkl"];
      }
    }
    if(!isArgumentGiven("PPL")){
      if (isArgumentGiven("CoulombFactors")) {

        // Read the factorsSliceSize.
        auto LambdaGR(getTensorArgument<complex>("CoulombFactors"));
        LambdaGR->set_name("LambdaGR");

        int NR(LambdaGR->lens[1]);

        int factorsSliceSize(
          getIntegerArgument("factorsSliceSize", DEFAULT_SLICE_SIZE)
        );
        if (factorsSliceSize == -1) {
          if (isArgumentGiven("factorsSliceFactor")) {
            double factorsSliceFactor(getRealArgument("factorsSliceFactor"));
            factorsSliceSize = NR * factorsSliceFactor;
          } else {
            factorsSliceSize = Nv;
          }
        }

        // Slice loop starts here
        for (int b(0); b < NR; b += factorsSliceSize) {
          for (int a(0); a < NR; a += factorsSliceSize) {
            LOG(1, getCapitalizedAbbreviation()) <<
              "Evaluting Fabij at R=" << a << ", S=" << b << std::endl;
            auto Fabij(
              sliceAmplitudesFromCoupledCoulombFactors(
                amplitudes, a, b, factorsSliceSize
              )
            );
            Fabij->set_name("Fabij");
            (*Rabij)["abij"] += (*Fabij)["abij"];
            delete Fabij;
          }
        }
      } else {
        // Read the integralsSliceSize. If not provided use No
        int integralsSliceSize(getIntegerArgument("integralsSliceSize",DEFAULT_SLICE_SIZE));
        if (integralsSliceSize == -1) {
          if (isArgumentGiven("integralsSliceFactor")) {
            double integralsSliceFactor(getRealArgument("integralsSliceFactor"));
            integralsSliceSize = Nv * integralsSliceFactor;
          } else {
            integralsSliceSize = No;
          }
        }

        Tensor<> realDressedGammaGab(realGammaGab);
        Tensor<> imagDressedGammaGab(imagGammaGab);
        realDressedGammaGab.set_name("realDressedGammaGab");
        imagDressedGammaGab.set_name("imagDressedGammaGab");
        // Construct dressed Coulomb vertex GammaGab
        realDressedGammaGab["Gab"] += (-1.0) * realGammaGai["Gbk"] * (*Tai)["ak"];
        imagDressedGammaGab["Gab"] += (-1.0) * imagGammaGai["Gbk"] * (*Tai)["ak"];

        // Slice loop starts here
        for (int b(0); b < Nv; b += integralsSliceSize) {
          // Slice the right GammaGab
 	  int rightGammaStart[] = { 0, b, 0};
          int rightGammaEnd[] = { NG, std::min(b+integralsSliceSize, Nv), Nv};
          auto realRightGamma(realDressedGammaGab.slice(rightGammaStart, rightGammaEnd));
          auto imagRightGamma(imagDressedGammaGab.slice(rightGammaStart, rightGammaEnd));

          for (int a(b); a < Nv; a += integralsSliceSize) {

            LOG(1, getCapitalizedAbbreviation()) <<
              "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
            // Slice the left GammaGab
            int leftGammaStart[] = { 0, a, 0};
            int leftGammaEnd[] = {NG, std::min(a+integralsSliceSize, Nv), Nv};
            auto realLeftGamma(realDressedGammaGab.slice(leftGammaStart, leftGammaEnd));
            auto imagLeftGamma(imagDressedGammaGab.slice(leftGammaStart, leftGammaEnd));

            int lenscd[] = {
              realLeftGamma.lens[1], realRightGamma.lens[1], realLeftGamma.lens[2], realRightGamma.lens[2]
            };
            int syms[] = {NS, NS, NS, NS};
            auto Vxycd(new Tensor<>(4, lenscd, syms, *realDressedGammaGab.wrld, "Vxycd"));
            (*Vxycd)["xycd"]  = realLeftGamma["Gxc"] * realRightGamma["Gyd"];
            (*Vxycd)["xycd"] += imagLeftGamma["Gxc"] * imagRightGamma["Gyd"];

            int lensij[] = { Vxycd->lens[0], Vxycd->lens[1], No, No};
            Tensor<> Rxyij(4, lensij, syms, *Vxycd->wrld, "Rxyij");

            // Contract sliced Vxycd with T2 and T1 Amplitudes using Xabij
            Rxyij["xyij"] = (*Vxycd)["xycd"] * Xabij["cdij"];

            sliceIntoResiduum(Rxyij, a, b, *Rabij);
            // The integrals of this slice are not needed anymore
            delete Vxycd;
          }
        }
      }
    }
    //********************************************************************************
    //***********************  T1 amplitude equations  *******************************
    //********************************************************************************
    {
      LOG(1, getCapitalizedAbbreviation()) <<
        "Solving T1 Amplitude Equations" << std::endl;

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
    }
  }
  return residuum;
}


PTR(FockVector<cc4s::complex>) CcsdEnergyFromCoulombIntegrals::getResiduum(
  const int i, const PTR(const FockVector<cc4s::complex>) &amplitudes
) {
  // Read Vabij integrals
  auto Vabij(getTensorArgument<complex>("PPHHCoulombIntegrals"));

  // get singles and doubles part of the amplitudes
  auto Tai( amplitudes->get(0) );
  Tai->set_name("Tai");
  auto Tabij( amplitudes->get(1) );
  Tabij->set_name("Tabij");

  // create residuum and get their singles and doubles part
  auto residuum( NEW(FockVector<complex>, *amplitudes) );
  *residuum *= complex(0.0);
  auto Rai( residuum->get(0) );
  Rai->set_name("Rai");
  auto Rabij( residuum->get(1) );
  Rabij->set_name("Rabij");

  if (i == 0 && !isArgumentGiven("initialDoublesAmplitudes") ) {
    // For first iteration compute only the MP2 amplitudes
    // Since Tabij = 0, Vabij is the only non-zero term
    LOG(1, getCapitalizedAbbreviation()) << "MP2 T2 Amplitudes" << std::endl;
    (*Rabij)["abij"] = (*Vabij)["abij"];
  } else {
    // For the rest iterations compute the CCSD amplitudes

    // Read all required integrals
    auto Vaijb(getTensorArgument<complex>("PHHPCoulombIntegrals"));
    auto Vijab(getTensorArgument<complex>("HHPPCoulombIntegrals"));
    auto Vaibj(getTensorArgument<complex>("PHPHCoulombIntegrals"));
    auto Vijkl(getTensorArgument<complex>("HHHHCoulombIntegrals"));
    auto Vijka(getTensorArgument<complex>("HHHPCoulombIntegrals"));
    auto Vaijk(getTensorArgument<complex>("PHHHCoulombIntegrals"));

    // Read the Coulomb vertex GammaGqr
    auto GammaGqr( getTensorArgument<complex>("CoulombVertex"));

    // Compute the No,Nv,NG,Np
    int No(Vabij->lens[2]);
    int Nv(Vabij->lens[0]);
    int NG(GammaGqr->lens[0]);
    int Np(GammaGqr->lens[1]);

    // Allocate and compute GammaGab,GammaGai,GammaGij from GammaGqr
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

    std::array<int,4> syms({{ NS, NS, NS, NS }});
    std::array<int,4> voov({{ Nv, No, No, Nv }});
    std::array<int,2> vv({{ Nv, Nv }});
    std::array<int,2> vo({{ Nv, No }});
    std::array<int,2> oo({{ No, No }});

    int distinguishable(
      getIntegerArgument("distinguishable", DEFAULT_DISTINGUISHABLE)
    );

    //********************************************************************************
    //***********************  T2 amplitude equations  *******************************
    //********************************************************************************

    LOG(1, getCapitalizedAbbreviation()) <<
      "Solving T2 Amplitude Equations" << std::endl;

    // Intermediates used both by T1 and T2
    Tensor<complex> Kac(2, vv.data(), syms.data(), *Vabij->wrld, "Kac");
    Tensor<complex> Kki(2, oo.data(), syms.data(), *Vabij->wrld, "Kki");
    Tensor<complex> Xabij(*Tabij);
    Xabij.set_name("Xabij");
    Xabij["abij"] += (*Tai)["ai"] * (*Tai)["bj"];

    {
      // Intermediates used for T2 amplitudes
      Tensor<complex> Lac(2, vv.data(), syms.data(), *Vabij->wrld, "Lac");
      Tensor<complex> Lki(2, oo.data(), syms.data(), *Vabij->wrld, "Lki");

      // Build Kac
      Kac["ac"]  = (-2.0) * (*Vijab)["klcd"] * Xabij["adkl"];
      Kac["ac"] += ( 1.0) * (*Vijab)["kldc"] * Xabij["adkl"];

      // Build Lac
      if (!distinguishable) {
        Lac["ac"]  = Kac["ac"];
      } else {
        // Intermediate tensor Yabij=T2-2*T1*T1
        Tensor<complex> Yabij(*Tabij);
        Yabij.set_name("Yabij");
        Yabij["abij"] += ( 2.0) * (*Tai)["ai"] * (*Tai)["bj"];
        Lac["ac"]      = (-1.0) * (*Vijab)["klcd"] * Yabij["adkl"]; // Use Yabij in DCSD
        Lac["ac"]     += ( 0.5) * (*Vijab)["kldc"] * Yabij["adkl"]; // Use Yabij in DCSD
      }
      Lac["ac"] += ( 2.0) * conjTransposeGammaGab["Gac"] * (*GammaGia)["Gkd"] * (*Tai)["dk"];
      Lac["ac"] += (-1.0) * conjTransposeGammaGab["Gad"] * (*GammaGia)["Gkc"] * (*Tai)["dk"];

      // Build Kki
      Kki["ki"]  = ( 2.0) * (*Vijab)["klcd"] * Xabij["cdil"];
      Kki["ki"] += (-1.0) * (*Vijab)["kldc"] * Xabij["cdil"];

      // Build Lki
      if (!distinguishable) {
        Lki["ki"]  = ( 1.0) *   Kki   ["ki"];
      } else {
        // Intermediate tensor Yabij=T2-2*T1*T1
        Tensor<complex> Yabij(*Tabij);
        Yabij.set_name("Yabij");
        Yabij["abij"] += ( 2.0) * (*Tai)["ai"] * (*Tai)["bj"];
        Lki["ki"]      = ( 1.0) * (*Vijab)["klcd"] * Yabij["cdil"]; // Use Yabij in DCSD
        Lki["ki"]     += (-0.5) * (*Vijab)["kldc"] * Yabij["cdil"]; // Use Yabij in DCSD
      }
      Lki["ki"] += ( 2.0) * (*Vijka)["klic"] * (*Tai)["cl"];
      Lki["ki"] += (-1.0) * (*Vijka)["lkic"] * (*Tai)["cl"];

      // Contract Lac with T2 Amplitudes
      (*Rabij)["abij"] += ( 1.0) * Lac["ac"] * (*Tabij)["cbij"];

      // Contract Lki with T2 Amplitudes
      (*Rabij)["abij"] += (-1.0) * Lki["ki"] * (*Tabij)["abkj"];
    }

    {
      // Contract Coulomb integrals with T2 amplitudes
      Tensor<complex> conjTransposeDressedGammaGai(conjTransposeGammaGai);
      conjTransposeDressedGammaGai.set_name("conjTransposeDressedGammaGai");
      conjTransposeDressedGammaGai["Gai"] += (-1.0) * conjTransposeGammaGij["Gki"] * (*Tai)["ak"];
      (*Rabij)["abij"] += ( 1.0) * conjTransposeDressedGammaGai["Gai"] * (*GammaGab)["Gbc"] * (*Tai)["cj"];

      (*Rabij)["abij"] += (-1.0) * (*Vaijk)["akij"] * (*Tai)["bk"];
      (*Rabij)["abij"] += (-1.0) * (*Tai)["bk"] * (*Vaijb)["akic"] * (*Tai)["cj"];
    }

    {
      // Build Xakic
      Tensor<complex> Xakic(4, voov.data(), syms.data(), *Vabij->wrld, "Xakic");

      Tensor<complex> conjTransposeDressedGammaGai(conjTransposeGammaGai);
      conjTransposeDressedGammaGai.set_name("conjTransposeDressedGammaGai");

      conjTransposeDressedGammaGai["Gai"] += (-1.0) * conjTransposeGammaGij["Gli"] * (*Tai)["al"];

      conjTransposeDressedGammaGai["Gai"] += ( 1.0) * conjTransposeGammaGab["Gad"] * (*Tai)["di"];

      // Intermediate tensor Yabij=T2-2*T1*T1
      Tensor<complex> Yabij(*Tabij);
      Yabij.set_name("Yabij");
      Yabij["abij"] += ( 2.0) * (*Tai)["ai"] * (*Tai)["bj"];

      conjTransposeDressedGammaGai["Gai"] += (-0.5) * conjTransposeGammaGia["Gld"] * Yabij["dail"];

      Xakic["akic"]  = ( 1.0) * conjTransposeDressedGammaGai["Gai"] * (*GammaGia)["Gkc"];

      Yabij["dclk"]  = ( 1.0) * (*Vijab)["lkdc"];
      if (!distinguishable) {
        Yabij["dclk"] += (-0.5) * (*Vijab)["lkcd"];
      }
      Xakic["akic"] += Yabij["dclk"] * (*Tabij)["adil"];

      // Contract and Xakic intermediates with T2 amplitudes Tabij
      Yabij["cbkj"]  = ( 2.0) * (*Tabij)["cbkj"];
      Yabij["cbkj"] += (-1.0) * (*Tabij)["bckj"];

      (*Rabij)["abij"] += ( 1.0) * Xakic["akic"] * Yabij["cbkj"];
    }

    {
      // Build Xakci
      Tensor<complex> Xakci(false, *Vaibj);
      Xakci.set_name("Xakci");

      // Construct dressed Coulomb vertex GammaGab and GammaGij
      Tensor<complex> conjTransposeDressedGammaGab(conjTransposeGammaGab);
      conjTransposeDressedGammaGab.set_name("conjTransposeDressedGammaGab");

      Tensor<complex> DressedGammaGij(*GammaGij);
      DressedGammaGij.set_name("DressedGammaGij");

      conjTransposeDressedGammaGab["Gac"] += (-1.0) * conjTransposeGammaGia["Glc"] * (*Tai)["al"];

      DressedGammaGij["Gki"] += ( 1.0) * (*GammaGia)["Gkd"] * (*Tai)["di"];

      // Xakci = Vakci - Vlkci * Tal + Vakcd * Tdi - Vlkcd * Tdail
      Xakci["akci"]  = ( 1.0) * conjTransposeDressedGammaGab["Gac"] * DressedGammaGij["Gki"];


      // Xakci = 0.5 * Vlkcd * Tdail
      if (!distinguishable) {
        Xakci["akci"] += (-0.5) * (*Vijab)["lkcd"] * (*Tabij)["dail"];
      }

      // Contract and Xakci intermediates with T2 amplitudes Tabij
      (*Rabij)["abij"] += (-1.0) * Xakci["akci"] * (*Tabij)["cbkj"];
      (*Rabij)["abij"] += (-1.0) * Xakci["bkci"] * (*Tabij)["ackj"];

      // Symmetrize Rabij by applying permutation operator
      (*Rabij)["abij"] += (*Rabij)["baji"];
    }

    //////////////////////////////////////////////////////////////////////
    // Now add all terms to Rabij that do not need to be symmetrized with
    // the permutation operator
    //////////////////////////////////////////////////////////////////////

    // Add Vabij to Rabij (MP2 term)
    (*Rabij)["abij"] += (*Vabij)["abij"];

    {
      // Build Xklij intermediate
      Tensor<complex> Xklij(false, *Vijkl);
      Xklij.set_name("Xklij");

      Xklij["klij"]  = (*Vijkl)["klij"];
      Xklij["klij"] += (*Vijka)["klic"] * (*Tai)["cj"];
      Xklij["klij"] += (*Vijka)["lkjc"] * (*Tai)["ci"];

      // Contract Xklij with T2+T1*T1 Amplitudes via Xabij
      (*Rabij)["abij"] +=  Xklij["klij"] * Xabij["abkl"];

      // Construct last term
      if (!distinguishable) {
        Xklij["klij"]  = (*Vijab)["klcd"] * Xabij["cdij"];
      } else {
        Xklij["klij"]  = (*Tai)["ci"] * (*Vijab)["klcd"] * (*Tai)["dj"];
      }

      // Add last term contracted only with the doubles
      // The singles term is computed in the slicing
      (*Rabij)["abij"] +=  Xklij["klij"] * (*Tabij)["abkl"];
    }

    if (isArgumentGiven("CoulombFactors")) {

      // Read the factorsSliceSize.
      auto LambdaGR(getTensorArgument<complex>("CoulombFactors"));
      LambdaGR->set_name("LambdaGR");

      int NR(LambdaGR->lens[1]);

      int factorsSliceSize(
        getIntegerArgument("factorsSliceSize", DEFAULT_SLICE_SIZE)
      );
      if (factorsSliceSize == -1) {
        if (isArgumentGiven("factorsSliceFactor")) {
          double factorsSliceFactor(getRealArgument("factorsSliceFactor"));
          factorsSliceSize = NR * factorsSliceFactor;
        } else {
          factorsSliceSize = Nv;
        }
      }

      // Slice loop starts here
      for (int b(0); b < NR; b += factorsSliceSize) {
        for (int a(0); a < NR; a += factorsSliceSize) {
          LOG(1, getCapitalizedAbbreviation()) <<
            "Evaluting Fabij at R=" << a << ", S=" << b << std::endl;
          auto Fabij(
            sliceAmplitudesFromCoupledCoulombFactors(
              amplitudes, a, b, factorsSliceSize
            )
          );
          Fabij->set_name("Fabij");
          (*Rabij)["abij"] += (*Fabij)["abij"];
          delete Fabij;
        }
      }
    } else {
      // Read the integralsSliceSize. If not provided use No
      int integralsSliceSize(getIntegerArgument("integralsSliceSize",DEFAULT_SLICE_SIZE));
      if (integralsSliceSize == -1) {
        if (isArgumentGiven("integralsSliceFactor")) {
          double integralsSliceFactor(getRealArgument("integralsSliceFactor"));
          integralsSliceSize = Nv * integralsSliceFactor;
        } else {
          integralsSliceSize = No;
        }
      }

      // Slice loop starts here
      for (int b(0); b < Nv; b += integralsSliceSize) {
        for (int a(b); a < Nv; a += integralsSliceSize) {
          LOG(1, getCapitalizedAbbreviation()) <<
            "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
          auto Vxycd(
            sliceCoupledCoulombIntegrals(amplitudes, a, b, integralsSliceSize)
          );
          Vxycd->set_name("Vxycd");
          int lens[] = { Vxycd->lens[0], Vxycd->lens[1], No, No };
          int syms[] = {NS, NS, NS, NS};
          Tensor<complex> Rxyij(4, lens, syms, *Vxycd->wrld, "Rxyij");

          // Contract sliced Vxycd with T2 and T1 Amplitudes using Xabij
          Rxyij["xyij"] = (*Vxycd)["xycd"] * Xabij["cdij"];

          sliceIntoResiduum(Rxyij, a, b, *Rabij);
          // The integrals of this slice are not needed anymore
          delete Vxycd;
        }
      }

    }

    //********************************************************************************
    //***********************  T1 amplitude equations  *******************************
    //********************************************************************************
    {
      LOG(1, getCapitalizedAbbreviation()) <<
        "Solving T1 Amplitude Equations" << std::endl;

      // Intermediates used for T1 amplitudes
      Tensor<complex> Kck(2, vo.data(), syms.data(), *Vabij->wrld, "Kck");

      // Contract Kac and Kki with T1 amplitudes
      (*Rai)["ai"] += ( 1.0) * Kac["ac"] * (*Tai)["ci"];
      (*Rai)["ai"] += (-1.0) * Kki["ki"] * (*Tai)["ak"];

      // Build Kck
      Kck["ck"]  = ( 2.0) * (*Vijab)["klcd"] * (*Tai)["dl"];
      Kck["ck"] += (-1.0) * (*Vijab)["kldc"] * (*Tai)["dl"];

      // Contract all the rest terms with T1 and T2 amplitudes
      (*Rai)["ai"] += ( 2.0) * Kck["ck"] * (*Tabij)["caki"];
      (*Rai)["ai"] += (-1.0) * Kck["ck"] * (*Tabij)["caik"];
      (*Rai)["ai"] += ( 1.0) * (*Tai)["ak"] * Kck["ck"] * (*Tai)["ci"];
      (*Rai)["ai"] += ( 2.0) * (*Vaijb)["akic"] * (*Tai)["ck"];
      (*Rai)["ai"] += (-1.0) * (*Vaibj)["akci"] * (*Tai)["ck"];

      (*Rai)["ai"] += ( 2.0) * conjTransposeGammaGab["Gac"] * (*GammaGia)["Gkd"] * Xabij["cdik"];
      (*Rai)["ai"] += (-1.0) * conjTransposeGammaGab["Gad"] * (*GammaGia)["Gkc"] * Xabij["cdik"];

      (*Rai)["ai"] += (-2.0) * (*Vijka)["klic"] * Xabij["ackl"];
      (*Rai)["ai"] += ( 1.0) * (*Vijka)["lkic"] * Xabij["ackl"];
    }

    delete GammaGij;
    delete GammaGia;
    delete GammaGai;
    delete GammaGab;
  }
  return residuum;
}

