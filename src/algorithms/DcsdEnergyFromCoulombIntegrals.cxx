#include <algorithms/DcsdEnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DcsdEnergyFromCoulombIntegrals);

DcsdEnergyFromCoulombIntegrals::DcsdEnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): ClusterSinglesDoublesAlgorithm(argumentList) {
}

DcsdEnergyFromCoulombIntegrals::~DcsdEnergyFromCoulombIntegrals() {
}

//////////////////////////////////////////////////////////////////////
// Hirata iteration routine for the CCSD amplitudes Tabij and Tai from
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
// modified to give DCSD amplitudes according to
// D. Kats, et. al., J. Chem. Phys. 142, 064111 (2015)
//////////////////////////////////////////////////////////////////////
void DcsdEnergyFromCoulombIntegrals::iterate(int i) {
  {
    // Read the amplitudes Tai and Tabij
    Tensor<> *Tai(&TaiMixer->getNext());
    Tai->set_name("Tai");
    Tensor<> *Tabij(&TabijMixer->getNext());
    Tabij->set_name("Tabij");

    // Read the Coulomb Integrals Vabcd Vabij Vaibj Vijkl Vabci Vijka
    // the PPPPCoulombIntegrals may not be given then slicing is required
    Tensor<> *Vabcd(isArgumentGiven("PPPPCoulombIntegrals") ?
		    getTensorArgument("PPPPCoulombIntegrals") : nullptr);
    Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
    Tensor<> *Vaibj(getTensorArgument("PHPHCoulombIntegrals"));
    Tensor<> *Vijkl(getTensorArgument("HHHHCoulombIntegrals"));
    Tensor<> *Vijka(getTensorArgument("HHHPCoulombIntegrals"));
    Tensor<> *Vabci(getTensorArgument("PPPHCoulombIntegrals"));

    // Read the Particle/Hole Eigenenergies epsi epsa
    Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
    Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));

    // Get abbreviation of algorithm
    std::string abbreviation(getAbbreviation());
    std::transform(abbreviation.begin(), abbreviation.end(), 
		   abbreviation.begin(), ::toupper);
  
    // Compute the No,Nv
    int No(epsi->lens[0]);
    int Nv(epsa->lens[0]);

    // Set symmetries for defining tensors
    int syms[] = { NS, NS, NS, NS };

    //********************************************************************************
    //***********************  T2 amplitude equations  *******************************
    //********************************************************************************

    {
    LOG(1, abbreviation) << "Solving T2 Amplitude Equations" << std::endl;

    // Allocate Tensors for T2 amplitudes
    Tensor<> Rabij(false, *Vabij);
    Rabij.set_name("Rabij");

    {
    // Intermediates used for T2 amplitudes
    int vv[] = { Nv, Nv };
    Tensor<> Lac(2, vv, syms, *epsi->wrld, "Lac");
    int oo[] = { No, No };
    Tensor<> Lki(2, oo, syms, *epsi->wrld, "Lki");

    Tensor<> Xklij(false, *Vijkl);
    Xklij.set_name("Xklij");
    Tensor<> Xakci(false, *Vaibj);
    Xakci.set_name("Xakci");
    int voov[] = { Nv, No, No, Nv };
    Tensor<> Xakic(4, voov, syms, *epsi->wrld, "Xakic");

    // Intermediates for quadratic T1 amplitudes contraction
    int vo[] = { Nv, No };
    Tensor<> Yai(2, vo, syms, *epsi->wrld, "Yai");
    Tensor<> Yijka(false, *Vijka);
    Yijka.set_name("Yijka");

    // Build Lac
    Lac["ac"]  = -1.0 * (*Vabij)["cdkl"] * (*Tabij)["adkl"]; // Multiplied by 0.5 in DCSD
    Lac["ac"] +=  0.5 * (*Vabij)["dckl"] * (*Tabij)["adkl"]; // Multiplied by 0.5 in DCSD
    Yai["ck"]  =  (*Vabij)["cdkl"] * (*Tai)["dl"];
    Lac["ac"] -=  2.0 * Yai["ck"] * (*Tai)["ak"];
    //Lac["ac"] -=  2.0 * (*Vabij)["cdkl"] * (*Tai)["ak"] * (*Tai)["dl"];
    Yai["ck"]  =  (*Vabij)["dckl"] * (*Tai)["dl"];
    Lac["ac"] +=  Yai["ck"] * (*Tai)["ak"];
    //Lac["ac"] +=  1.0 * (*Vabij)["dckl"] * (*Tai)["ak"] * (*Tai)["dl"];
    Lac["ac"] +=  2.0 * (*Vabci)["cdak"] * (*Tai)["dk"];
    Lac["ac"] -=  (*Vabci)["dcak"] * (*Tai)["dk"];

    // Build Lki
    Lki["ki"]  = 1.0 * (*Vabij)["cdkl"] * (*Tabij)["cdil"]; // Multiplied by 0.5 in DCSD
    Lki["ki"] -= 0.5 * (*Vabij)["dckl"] * (*Tabij)["cdil"]; // Multiplied by 0.5 in DCSD
    Yai["ck"]  = (*Vabij)["cdkl"] * (*Tai)["dl"];
    Lki["ki"] += 2.0 * Yai["ck"] * (*Tai)["ci"];
    //Lki["ki"] += 2.0 * (*Vabij)["cdkl"] * (*Tai)["ci"] * (*Tai)["dl"];
    Yai["ck"]  = (*Vabij)["dckl"] * (*Tai)["dl"];
    Lki["ki"] -= Yai["ck"] * (*Tai)["ci"];
    //Lki["ki"] -= 1.0 * (*Vabij)["dckl"] * (*Tai)["ci"] * (*Tai)["dl"];
    Lki["ki"] += 2.0 * (*Vijka)["klic"] * (*Tai)["cl"];
    Lki["ki"] -= (*Vijka)["lkic"] * (*Tai)["cl"];
    
    // Contract Lac with T2 Amplitudes
    Rabij["abij"] = Lac["ac"] * (*Tabij)["cbij"];

    // Contract Lki with T2 Amplitudes
    Rabij["abij"] -= Lki["ki"] * (*Tabij)["abkj"];

    // Contract Coulomb integrals with T2 amplitudes
    Rabij["abij"] += (*Vabci)["baci"] * (*Tai)["cj"];
    Yijka["kjib"]  = (*Vaibj)["bkci"] * (*Tai)["cj"];
    Rabij["abij"] -= Yijka["kjib"] * (*Tai)["ak"];
    //Rabij["abij"] -= (*Vaibj)["bkci"] * (*Tai)["ak"] * (*Tai)["cj"];
    Rabij["abij"] -= (*Vijka)["jika"] * (*Tai)["bk"];
    Yijka["jika"]  = (*Vabij)["acik"] * (*Tai)["cj"];
    Rabij["abij"] += Yijka["jika"] * (*Tai)["bk"];
    //Rabij["abij"] += (*Vabij)["acik"] * (*Tai)["cj"] * (*Tai)["bk"];

    // Build Xakic
    Xakic["akic"]  = (*Vabij)["acik"];
    Xakic["akic"] -= (*Vijka)["lkic"] * (*Tai)["al"];
    Xakic["akic"] += (*Vabci)["acdk"] * (*Tai)["di"];
    Xakic["akic"] -= 0.5 * (*Vabij)["dclk"] * (*Tabij)["dail"];
    Yijka["iklc"]  = (*Vabij)["dclk"] * (*Tai)["di"];
    Xakic["akic"] -= Yijka["iklc"] * (*Tai)["al"];
    //Xakic["akic"] -= (*Vabij)["dclk"] * (*Tai)["di"] * (*Tai)["al"];
    Xakic["akic"] += (*Vabij)["dclk"] * (*Tabij)["adil"];
    //Xakic["akic"] -= 0.5 * (*Vabij)["cdlk"] * (*Tabij)["adil"]; // Removed in DCSD

    // Build Xakci
    Xakci["akci"]  = (*Vaibj)["akci"];
    Xakci["akci"] -= (*Vijka)["klic"] * (*Tai)["al"];
    Xakci["akci"] += (*Vabci)["adck"] * (*Tai)["di"];
    //Xakci["akci"] -= 0.5 * (*Vabij)["cdlk"] * (*Tabij)["dail"]; // Removed in DCSD
    Yijka["ilkc"]  = (*Vabij)["cdlk"] * (*Tai)["di"];
    Xakci["akci"] -= Yijka["ilkc"] * (*Tai)["al"];
    //Xakci["akci"] -= (*Vabij)["cdlk"] * (*Tai)["di"] * (*Tai)["al"];

    // Contract Xakic and Xakci intermediates with T2 amplitudes Tabij
    Rabij["abij"] += 2.0 * Xakic["akic"] * (*Tabij)["cbkj"];
    Rabij["abij"] -= Xakic["akic"] * (*Tabij)["bckj"];

    Rabij["abij"] -= Xakci["akci"] * (*Tabij)["cbkj"];
    Rabij["abij"] -= Xakci["bkci"] * (*Tabij)["ackj"];

    // Symmetrize Rabij by applying permutation operator
    // to save memory we use Xakci as intermediate for the permutation operator 
    Xakci["aibj"]  = Rabij["abij"];
    Rabij["abij"] += Xakci["bjai"]; 

    //////////////////////////////////////////////////////////////////////
    // Now add all terms to Rabij that do not need to be symmetrized with
    // the permutation operator
    //////////////////////////////////////////////////////////////////////

    // Rabij are the Tabij amplitudes for the next iteration and need to be build
    Rabij["abij"] += (*Vabij)["abij"];

    // Build Xklij intermediate
    Xklij["klij"]  = (*Vijkl)["klij"];
    Xklij["klij"] += (*Vijka)["klic"] * (*Tai)["cj"];
    Xklij["klij"] += (*Vijka)["lkjc"] * (*Tai)["ci"];
    Yijka["jklc"]  = (*Vabij)["cdkl"] * (*Tai)["dj"]; 
    Xklij["klij"] += Yijka["jklc"] * (*Tai)["ci"];
    //Xklij["klij"] += (*Vabij)["cdkl"] * (*Tai)["ci"] * (*Tai)["dj"]; 

    // Contract Xklij with T2 Amplitudes
    Rabij["abij"] += Xklij["klij"] * (*Tabij)["abkl"];

    // Contract Xklij with T1 Amplitudes
    Xklij["klij"] += (*Vabij)["cdkl"] * (*Tabij)["cdij"]; //Removed in Dcsd from T2 Amplitudes
    Yijka["lija"]  = Xklij["klij"] * (*Tai)["ak"];
    Rabij["abij"] += Yijka["lija"] * (*Tai)["bl"];
    //Rabij["abij"] += Xklij["klij"] * (*Tai)["ak"] * (*Tai)["bl"];
    }

    if (Vabcd) {
      // Build Xabcd intermediate
      Tensor<> Xabcd(Vabcd);
      Xabcd.set_name("Xabcd");
      Xabcd["abcd"] -= (*Vabci)["cdak"] * (*Tai)["bk"];
      Xabcd["abcd"] -= (*Vabci)["dcbk"] * (*Tai)["ak"];

      // Construct intermediate tensor
      Tensor<> Xabij(Tabij);
      Xabij.set_name("Xabij");
      Xabij["abij"] += (*Tai)["ai"] * (*Tai)["bj"];

      // Contract Xabcd with T2 and T1 Amplitudes using Xabij
      Rabij["abij"] += Xabcd["abcd"] * Xabij["cdij"];

      /*
      // Contract Xabcd with T2 and T1 Amplitudes
      Rabij["abij"] += Xabcd["abcd"] * (*Tabij)["cdij"];
      Rabij["abij"] += Xabcd["abcd"] * (*Tai)["ci"] * (*Tai)["dj"];
      */
    } else {
      // Slice if Vabcd is not specified

      // Read the sliceRank. If not provided use No
      int sliceRank(getIntegerArgument
		    ("sliceRank",No));

      // Slice loop starts here
      for (int b(0); b < Nv; b += sliceRank) {
        for (int a(b); a < Nv; a += sliceRank) {
          LOG(1, abbreviation) << "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
          // get the sliced integrals already coupled to the singles
          Tensor<> *Xxycd(sliceCoupledCoulombIntegrals(a, b, sliceRank));
	  Xxycd->set_name("Xxycd");
          int lens[] = { Xxycd->lens[0], Xxycd->lens[1], No, No };
          int syms[] = {NS, NS, NS, NS};
          Tensor<> Rxyij(4, lens, syms, *Xxycd->wrld, "Rxyij");
          Rxyij["xyij"] =  (*Xxycd)["xycd"] * (*Tabij)["cdij"];
          Rxyij["xyij"] += (*Xxycd)["xycd"] * (*Tai)["ci"] * (*Tai)["dj"];
          sliceIntoResiduum(Rxyij, a, b, Rabij);
          // the integrals of this slice are not needed anymore
          delete Xxycd;
        }
      }
    }

    // calculate the amplitdues from the residuum
    doublesAmplitudesFromResiduum(Rabij);
    // and append them to the mixer
    TabijMixer->append(Rabij);
    }

    //********************************************************************************
    //***********************  T1 amplitude equations  *******************************
    //********************************************************************************
    {
    LOG(1, abbreviation) << "Solving T1 Amplitude Equations" << std::endl;

    // Allocate Tensors for T1 amplitudes
    Tensor<> Rai(false, *Tai);
    Rai.set_name("Rai");

    // Intermediates used for T1 amplitudes
    int vo[] = { Nv, No };
    Tensor<> Kck(2, vo, syms, *epsi->wrld, "Kck");

    // Intermediates used for T1 amplitudes
    int vv[] = { Nv, Nv };
    Tensor<> Kac(2, vv, syms, *epsi->wrld, "Kac");
    int oo[] = { No, No };
    Tensor<> Kki(2, oo, syms, *epsi->wrld, "Kki");

    // Intermediate for quadratic T1 contractions
    Tensor<> Yij(2, oo, syms, *epsi->wrld, "Yij");
    Tensor<> Yai(2, vo, syms, *epsi->wrld, "Yai");

    // Build Kac
    Kac["ac"]  = -2.0 * (*Vabij)["cdkl"] * (*Tabij)["adkl"];
    Kac["ac"] += (*Vabij)["dckl"] * (*Tabij)["adkl"];
    Yai["ck"]  = (*Vabij)["cdkl"] * (*Tai)["dl"];
    Kac["ac"] -= 2.0 * Yai["ck"] * (*Tai)["ak"];
    //Kac["ac"] -= 2.0 * (*Vabij)["cdkl"] * (*Tai)["ak"] * (*Tai)["dl"];
    Yai["ck"]  = (*Vabij)["dckl"] * (*Tai)["dl"];
    Kac["ac"] += Yai["ck"] * (*Tai)["ak"];
    //Kac["ac"] += (*Vabij)["dckl"] * (*Tai)["ak"] * (*Tai)["dl"];

    // Build Kki
    Kki["ki"]  = 2.0 * (*Vabij)["cdkl"] * (*Tabij)["cdil"];
    Kki["ki"] -= (*Vabij)["dckl"] * (*Tabij)["cdil"];
    Yai["ck"]  = (*Vabij)["cdkl"] * (*Tai)["dl"];
    Kki["ki"] += 2.0 * Yai["ck"] * (*Tai)["ci"];
    //Kki["ki"] += 2.0 * (*Vabij)["cdkl"] * (*Tai)["ci"] * (*Tai)["dl"];
    Yai["ck"]  = (*Vabij)["dckl"] * (*Tai)["dl"];
    Kki["ki"] -= Yai["ck"] * (*Tai)["ci"];
    //Kki["ki"] -= (*Vabij)["dckl"] * (*Tai)["ci"] * (*Tai)["dl"];

    // Contract Kac and Kki with T1 amplitudes
    Rai["ai"]  = Kac["ac"] * (*Tai)["ci"];
    Rai["ai"] -= Kki["ki"] * (*Tai)["ak"];

    // Build Kck
    Kck["ck"]  = 2.0 * (*Vabij)["cdkl"] * (*Tai)["dl"];
    Kck["ck"] -= (*Vabij)["cdlk"] * (*Tai)["dl"];

    // Contract all the rest terms with T1 and T2 amplitudes
    Rai["ai"] += 2.0 * Kck["ck"] * (*Tabij)["caki"];
    Rai["ai"] -= Kck["ck"] * (*Tabij)["caik"];
    Yij["ik"]  = Kck["ck"] * (*Tai)["ci"];
    Rai["ai"] += Yij["ik"] * (*Tai)["ak"];
    //Rai["ai"] += Kck["ck"] * (*Tai)["ci"] * (*Tai)["ak"];
    Rai["ai"] += 2.0 * (*Vabij)["acik"] * (*Tai)["ck"];
    Rai["ai"] -= (*Vaibj)["akci"] * (*Tai)["ck"];
    Rai["ai"] += 2.0 * (*Vabci)["cdak"] * (*Tabij)["cdik"];
    Rai["ai"] -= (*Vabci)["dcak"] * (*Tabij)["cdik"];
    Kac["ac"]  = (*Vabci)["cdak"] * (*Tai)["dk"];     // use Kac to save memory
    Rai["ai"] += 2.0 * Kac["ac"] * (*Tai)["ci"];      // use Kac to save memory
    //Rai["ai"] += 2.0 * (*Vabci)["cdak"] * (*Tai)["ci"] * (*Tai)["dk"];
    Kac["ac"]  = (*Vabci)["dcak"] * (*Tai)["dk"];     // use Kac to save memory
    Rai["ai"] -= Kac["ac"] * (*Tai)["ci"];            // use Kac to save memory
    //Rai["ai"] -= (*Vabci)["dcak"] * (*Tai)["ci"] * (*Tai)["dk"];
    Rai["ai"] -= 2.0 * (*Vijka)["klic"] * (*Tabij)["ackl"];
    Rai["ai"] += (*Vijka)["lkic"] * (*Tabij)["ackl"];
    Yij["ki"]  = (*Vijka)["klic"] * (*Tai)["cl"];
    Rai["ai"] -= 2.0 * Yij["ki"] * (*Tai)["ak"];
    //Rai["ai"] -= 2.0 * (*Vijka)["klic"] * (*Tai)["ak"] * (*Tai)["cl"];
    Yij["ki"]  = (*Vijka)["lkic"] * (*Tai)["cl"];
    Rai["ai"] += Yij["ki"] * (*Tai)["ak"];
    //Rai["ai"] += (*Vijka)["lkic"] * (*Tai)["ak"] * (*Tai)["cl"];

    singlesAmplitudesFromResiduum(Rai);
    TaiMixer->append(Rai);
    }
  }
}

