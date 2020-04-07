#include <algorithms/CcsdEnergyFromCoulombIntegralsReference.hpp>
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

ALGORITHM_REGISTRAR_DEFINITION(CcsdEnergyFromCoulombIntegralsReference);

CcsdEnergyFromCoulombIntegralsReference::CcsdEnergyFromCoulombIntegralsReference(
  std::vector<Argument> const &argumentList
): ClusterSinglesDoublesAlgorithm(argumentList) {
}

CcsdEnergyFromCoulombIntegralsReference::~CcsdEnergyFromCoulombIntegralsReference() {
}

//////////////////////////////////////////////////////////////////////
// Hirata iteration routine for the CCSD amplitudes Tabij and Tai from
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
//////////////////////////////////////////////////////////////////////

PTR(FockVector<double>) CcsdEnergyFromCoulombIntegralsReference::getResiduum(
  const int i, const PTR(const FockVector<double>) &amplitudes
) {
  // get singles and doubles part of the amplitudes
  auto Tai( amplitudes->get(0) );
  Tai->set_name("Tai");
  auto Tabij( amplitudes->get(1) );
  Tabij->set_name("Tabij");
  const bool onlyPPL(getIntegerArgument("onlyPPL", 0) == 1);

  // create residuum and get their singles and doubles part
  auto residuum( NEW(FockVector<double>, *amplitudes) );
  *residuum *= 0.0;
  auto Rai( residuum->get(0) );
  Rai->set_name("Rai");
  auto Rabij( residuum->get(1) );
  Rabij->set_name("Rabij");

  // get part of Coulomb integrals used whether the amplitudes are zero or not
  auto Vabij(getTensorArgument("PPHHCoulombIntegrals"));

  if (i == 0 && !isArgumentGiven("initialDoublesAmplitudes") && !onlyPPL)  {
    // For first iteration compute only the MP2 amplitudes
    // Since Tabij = 0, Vabij is the only non-zero term
    LOG(1, getCapitalizedAbbreviation()) << "MP2 T2 Amplitudes" << std::endl;
    (*Rabij)["abij"] = (*Vabij)["abij"];
  } else if (onlyPPL) {

    LOG(1, getCapitalizedAbbreviation())
        << "Considering only PPL diagrams" << std::endl;

    auto Vabcd(getTensorArgument("PPPPCoulombIntegrals"));
    auto Vabci(getTensorArgument("PPPHCoulombIntegrals"));
    Tensor<> Xabcd(false, *Vabcd);

    // Build Xabcd intermediate
    Xabcd["abcd"]  = ( 1.0) * (*Vabcd)["abcd"];
    Xabcd["abcd"] += (-1.0) * (*Vabci)["cdak"] * (*Tai)["bk"];
    Xabcd["abcd"] += (-1.0) * (*Vabci)["dcbk"] * (*Tai)["ak"];

    (*Rabij)["abij"]  = Xabcd["abcd"] * (*Tabij)["cdij"];
    (*Rabij)["abij"] += Xabcd["abcd"] * (*Tai)["ci"] * (*Tai)["dj"];

  } else {
    // For the rest iterations compute the CCSD amplitudes

    // Read all required integrals
    auto Vabcd(getTensorArgument("PPPPCoulombIntegrals"));
    auto Vaibj(getTensorArgument("PHPHCoulombIntegrals"));
    auto Vijkl(getTensorArgument("HHHHCoulombIntegrals"));
    auto Vijka(getTensorArgument("HHHPCoulombIntegrals"));
    auto Vabci(getTensorArgument("PPPHCoulombIntegrals"));

    // Compute the No,Nv,NG,Np
    int No(Vabij->lens[2]);
    int Nv(Vabij->lens[0]);

    std::array<int,4> syms({{ NS, NS, NS, NS }});
    std::array<int,4> voov({{ Nv, No, No, Nv }});
    std::array<int,2> vv({{ Nv, Nv }});
    std::array<int,2> vo({{ Nv, No }});
    std::array<int,2> oo({{ No, No }});

    //*************************************************************************
    //****************  T2 amplitude equations  *******************************
    //*************************************************************************

    LOG(1, getCapitalizedAbbreviation()) <<
      "Solving T2 Amplitude Equations" << std::endl;

    // Define intermediates
    Tensor<> Lac(2, vv.data(), syms.data(), *Vabij->wrld, "Lac");
    Tensor<> Kac(2, vv.data(), syms.data(), *Vabij->wrld, "Kac");
    Tensor<> Lki(2, oo.data(), syms.data(), *Vabij->wrld, "Lki");
    Tensor<> Kki(2, oo.data(), syms.data(), *Vabij->wrld, "Kki");
    Tensor<> Kck(2, vo.data(), syms.data(), *Vabij->wrld, "Kck");

    Tensor<> Xklij(false, *Vijkl);
    Tensor<> Xakci(false, *Vaibj);
    Tensor<> Xabcd(false, *Vabcd);
    Tensor<> Xakic(4, voov.data(), syms.data(), *Vabij->wrld, "Xakic");

    // Build Kac
    Kac["ac"]  = (-2.0) * (*Vabij)["cdkl"] * (*Tabij)["adkl"];
    Kac["ac"] += ( 1.0) * (*Vabij)["dckl"] * (*Tabij)["adkl"];
    Kac["ac"] += (-2.0) * (*Vabij)["cdkl"] * (*Tai)["ak"] * (*Tai)["dl"];
    Kac["ac"] += ( 1.0) * (*Vabij)["dckl"] * (*Tai)["ak"] * (*Tai)["dl"];

    // Build Lac
    Lac["ac"]  = Kac["ac"];
    Lac["ac"] += ( 2.0) * (*Vabci)["cdak"] * (*Tai)["dk"];
    Lac["ac"] += (-1.0) * (*Vabci)["dcak"] * (*Tai)["dk"];

    // Build Kki
    Kki["ki"]  = ( 2.0) * (*Vabij)["cdkl"] * (*Tabij)["cdil"];
    Kki["ki"] += (-1.0) * (*Vabij)["dckl"] * (*Tabij)["cdil"];
    Kki["ki"] += ( 2.0) * (*Vabij)["cdkl"] * (*Tai)["ci"] * (*Tai)["dl"];
    Kki["ki"] += (-1.0) * (*Vabij)["dckl"] * (*Tai)["ci"] * (*Tai)["dl"];

    // Build Lki
    Lki["ki"]  = Kki["ki"];
    Lki["ki"] += ( 2.0) * (*Vijka)["klic"] * (*Tai)["cl"];
    Lki["ki"] += (-1.0) * (*Vijka)["lkic"] * (*Tai)["cl"];
    
    // Contract Lac with T2 Amplitudes
    (*Rabij)["abij"] += ( 1.0) * Lac["ac"] * (*Tabij)["cbij"];

    // Contract Lki with T2 Amplitudes
    (*Rabij)["abij"] += (-1.0) * Lki["ki"] * (*Tabij)["abkj"];

    // Contract Coulomb integrals with T2 amplitudes
    (*Rabij)["abij"] += ( 1.0) * (*Vabci)["baci"] * (*Tai)["cj"];
    (*Rabij)["abij"] += (-1.0) * (*Vaibj)["bkci"] * (*Tai)["ak"] * (*Tai)["cj"];
    (*Rabij)["abij"] += (-1.0) * (*Vijka)["jika"] * (*Tai)["bk"];
    (*Rabij)["abij"] += (-1.0) * (*Vabij)["acik"] * (*Tai)["cj"] * (*Tai)["bk"];

    // Build Xakic
    Xakic["akic"]  = (*Vabij)["acik"];
    Xakic["akic"] += (-1.0) * (*Vijka)["lkic"] * (*Tai)["al"];
    Xakic["akic"] += ( 1.0) * (*Vabci)["acdk"] * (*Tai)["di"];
    Xakic["akic"] += (-0.5) * (*Vabij)["dclk"] * (*Tabij)["dail"];
    Xakic["akic"] += (-1.0) * (*Vabij)["dclk"] * (*Tai)["di"] * (*Tai)["al"];
    Xakic["akic"] += ( 1.0) * (*Vabij)["dclk"] * (*Tabij)["adil"];
    Xakic["akic"] += (-0.5) * (*Vabij)["cdlk"] * (*Tabij)["adil"];

    // Build Xakci
    Xakci["akci"]  = (*Vaibj)["akci"];
    Xakci["akci"] += (-1.0) * (*Vijka)["klic"] * (*Tai)["al"];
    Xakci["akci"] += ( 1.0) * (*Vabci)["adck"] * (*Tai)["di"];
    Xakci["akci"] += (-0.5) * (*Vabij)["cdlk"] * (*Tabij)["dail"];
    Xakci["akci"] += (-1.0) * (*Vabij)["cdlk"] * (*Tai)["di"] * (*Tai)["al"];

    // Contract Xakic and Xakci intermediates with T2 amplitudes Tabij
    (*Rabij)["abij"] += ( 2.0) * Xakic["akic"] * (*Tabij)["cbkj"];
    (*Rabij)["abij"] += (-1.0) * Xakic["akic"] * (*Tabij)["bckj"];

    (*Rabij)["abij"] += (-1.0) * Xakci["akci"] * (*Tabij)["cbkj"];
    (*Rabij)["abij"] += (-1.0) * Xakci["bkci"] * (*Tabij)["ackj"];

    // Symmetrize Rabij by applying permutation operator
    // to save memory we use Xakci as intermediate for the permutation operator 
    Xakci["aibj"]  = (*Rabij)["abij"];
    (*Rabij)["abij"] += Xakci["bjai"]; 

    //////////////////////////////////////////////////////////////////////
    // Now add all terms to Rabij that do not need to be symmetrized with
    // the permutation operator
    //////////////////////////////////////////////////////////////////////

    // Rabij are the Tabij amplitudes for the next iteration and need to be build
    (*Rabij)["abij"] += (*Vabij)["abij"];

    // Build Xklij intermediate
    Xklij["klij"]  = (*Vijkl)["klij"];
    Xklij["klij"] += (*Vijka)["klic"] * (*Tai)["cj"];
    Xklij["klij"] += (*Vijka)["lkjc"] * (*Tai)["ci"];
    Xklij["klij"] += (*Vabij)["cdkl"] * (*Tabij)["cdij"];
    Xklij["klij"] += (*Vabij)["cdkl"] * (*Tai)["ci"] * (*Tai)["dj"]; 

    // Contract Xklij with T2 Amplitudes
    (*Rabij)["abij"] += Xklij["klij"] * (*Tabij)["abkl"];

    // Contract Xklij with T1 Amplitudes
    (*Rabij)["abij"] += Xklij["klij"] * (*Tai)["ak"] * (*Tai)["bl"];

    // Build Xabcd intermediate
    Xabcd["abcd"]  = ( 1.0) * (*Vabcd)["abcd"];
    Xabcd["abcd"] += (-1.0) * (*Vabci)["cdak"] * (*Tai)["bk"];
    Xabcd["abcd"] += (-1.0) * (*Vabci)["dcbk"] * (*Tai)["ak"];
    
    // Contract Xabcd with T2 and T1 Amplitudes
    (*Rabij)["abij"] += Xabcd["abcd"] * (*Tabij)["cdij"];
    (*Rabij)["abij"] += Xabcd["abcd"] * (*Tai)["ci"] * (*Tai)["dj"];

    //********************************************************************************
    //***********************  T1 amplitude equations  *******************************
    //********************************************************************************
      LOG(1, getCapitalizedAbbreviation()) <<
        "Solving T1 Amplitude Equations" << std::endl;

    // Contract Kac and Kki with T1 amplitudes
    (*Rai)["ai"] += ( 1.0) * Kac["ac"] * (*Tai)["ci"];
    (*Rai)["ai"] += (-1.0) * Kki["ki"] * (*Tai)["ak"];

    //Build Kck
    Kck["ck"]  = ( 2.0) * (*Vabij)["cdkl"] * (*Tai)["dl"];
    Kck["ck"] += (-1.0) * (*Vabij)["cdlk"] * (*Tai)["dl"];

    // Contract all the rest terms with T1 and T2 amplitudes
    (*Rai)["ai"] += ( 2.0) * Kck["ck"] * (*Tabij)["caki"];
    (*Rai)["ai"] += (-1.0) * Kck["ck"] * (*Tabij)["caik"];
    // TODO: check whether it's +
    (*Rai)["ai"] += ( 1.0) * Kck["ck"] * (*Tai)["ci"] * (*Tai)["ak"];
    (*Rai)["ai"] += ( 2.0) * (*Vabij)["acik"] * (*Tai)["ck"];
    (*Rai)["ai"] += (-1.0) * (*Vaibj)["akci"] * (*Tai)["ck"];
    (*Rai)["ai"] += ( 2.0) * (*Vabci)["cdak"] * (*Tabij)["cdik"];
    (*Rai)["ai"] += (-1.0) * (*Vabci)["dcak"] * (*Tabij)["cdik"];
    (*Rai)["ai"] += ( 2.0) * (*Vabci)["cdak"] * (*Tai)["ci"] * (*Tai)["dk"];
    (*Rai)["ai"] += (-1.0) * (*Vabci)["dcak"] * (*Tai)["ci"] * (*Tai)["dk"];
    (*Rai)["ai"] += (-2.0) * (*Vijka)["klic"] * (*Tabij)["ackl"];
    (*Rai)["ai"] += ( 1.0) * (*Vijka)["lkic"] * (*Tabij)["ackl"];
    (*Rai)["ai"] += (-2.0) * (*Vijka)["klic"] * (*Tai)["ak"] * (*Tai)["cl"];
    (*Rai)["ai"] += ( 1.0) * (*Vijka)["lkic"] * (*Tai)["ak"] * (*Tai)["cl"];
  }
  return residuum;
}

PTR(FockVector<cc4s::complex>) CcsdEnergyFromCoulombIntegralsReference::getResiduum(
  const int i, const PTR(const FockVector<cc4s::complex>) &amplitudes
) {
  throw new EXCEPTION("This is not implemented");
}
