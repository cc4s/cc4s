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
    // Read the DCSD amplitudes Tai and Tabij
    Tensor<> *Tai(&TaiMixer->getNext());
    Tensor<> *Tabij(&TabijMixer->getNext());

    // Read the Coulomb Integrals Vabcd Vabij Vaibj Vijkl Vabci Vijka
    // the PPPPCoulombIntegrals may not be given then slicing is required
    Tensor<> *Vabcd(
      isArgumentGiven("PPPPCoulombIntegrals") ?
        getTensorArgument("PPPPCoulombIntegrals") : nullptr
    );
    Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
    Tensor<> *Vaibj(getTensorArgument("PHPHCoulombIntegrals"));
    Tensor<> *Vijkl(getTensorArgument("HHHHCoulombIntegrals"));
    Tensor<> *Vijka(getTensorArgument("HHHPCoulombIntegrals"));
    Tensor<> *Vabci(getTensorArgument("PPPHCoulombIntegrals"));

    // Read the Particle/Hole Eigenenergies epsi epsa
    Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
    Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  
    // Compute the No,Nv
    int No(epsi->lens[0]);
    int Nv(epsa->lens[0]);

    int syms[] = { NS, NS, NS, NS };
    int voov[] = { Nv, No, No, Nv };
    int vo[] = { Nv, No };
    int vv[] = { Nv, Nv };
    int oo[] = { No, No };

    // Allocate Tensors for T1 amplitudes
    Tensor<> Rai(false, *Tai);
    //Tensor<> Dai(false, *Tai);

    // Allocate Tensors for T2 amplitudes
    Tensor<> Rabij(false, *Vabij);

    // Define intermediates
    Tensor<> Lac(2, vv, syms, *epsi->wrld, "Lac");
    Tensor<> Kac(2, vv, syms, *epsi->wrld, "Kac");
    Tensor<> Lki(2, oo, syms, *epsi->wrld, "Lki");
    Tensor<> Kki(2, oo, syms, *epsi->wrld, "Kki");
    Tensor<> Kck(2, vo, syms, *epsi->wrld, "Kck");

    Tensor<> Xklij(false, *Vijkl);
    Tensor<> Xakci(false, *Vaibj);
    Tensor<> Xakic(4, voov, syms, *epsi->wrld, "Xakic");

    //********************************************************************************
    //***********************  T2 amplitude equations  *******************************
    //********************************************************************************

    LOG(1, "DCSD") << "Solving T2 DCSD Amplitude Equations  ...";

    // Build Kac
    Kac["ac"]  = -2.0 * (*Vabij)["cdkl"] * (*Tabij)["adkl"];
    Kac["ac"] += (*Vabij)["dckl"] * (*Tabij)["adkl"];
    Kac["ac"] -= 2.0 * (*Vabij)["cdkl"] * (*Tai)["ak"] * (*Tai)["dl"];
    Kac["ac"] += (*Vabij)["dckl"] * (*Tai)["ak"] * (*Tai)["dl"];

    // Build Lac
    Lac["ac"]  = -1.0 * (*Vabij)["cdkl"] * (*Tabij)["adkl"]; // Multiplied by 0.5 in DCSD
    Lac["ac"] +=  0.5 * (*Vabij)["dckl"] * (*Tabij)["adkl"]; // Multiplied by 0.5 in DCSD
    Lac["ac"] -=  2.0 * (*Vabij)["cdkl"] * (*Tai)["ak"] * (*Tai)["dl"];
    Lac["ac"] +=  1.0 * (*Vabij)["dckl"] * (*Tai)["ak"] * (*Tai)["dl"];
    Lac["ac"] +=  2.0 * (*Vabci)["cdak"] * (*Tai)["dk"];
    Lac["ac"] -=  (*Vabci)["dcak"] * (*Tai)["dk"];

    // Build Kki
    Kki["ki"]  = 2.0 * (*Vabij)["cdkl"] * (*Tabij)["cdil"];
    Kki["ki"] -= (*Vabij)["dckl"] * (*Tabij)["cdil"];
    Kki["ki"] += 2.0 * (*Vabij)["cdkl"] * (*Tai)["ci"] * (*Tai)["dl"];
    Kki["ki"] -= (*Vabij)["dckl"] * (*Tai)["ci"] * (*Tai)["dl"];

    // Build Lki
    Lki["ki"]  = 1.0 * (*Vabij)["cdkl"] * (*Tabij)["cdil"]; // Multiplied by 0.5 in DCSD
    Lki["ki"] -= 0.5 * (*Vabij)["dckl"] * (*Tabij)["cdil"]; // Multiplied by 0.5 in DCSD
    Lki["ki"] += 2.0 * (*Vabij)["cdkl"] * (*Tai)["ci"] * (*Tai)["dl"];
    Lki["ki"] -= 1.0 * (*Vabij)["dckl"] * (*Tai)["ci"] * (*Tai)["dl"];
    Lki["ki"] += 2.0 * (*Vijka)["klic"] * (*Tai)["cl"];
    Lki["ki"] -= (*Vijka)["lkic"] * (*Tai)["cl"];
    
    // Contract Lac with T2 Amplitudes
    Rabij["abij"] = Lac["ac"] * (*Tabij)["cbij"];

    // Contract Lki with T2 Amplitudes
    Rabij["abij"] -= Lki["ki"] * (*Tabij)["abkj"];

    // Contract Coulomb integrals with T2 amplitudes
    Rabij["abij"] += (*Vabci)["baci"] * (*Tai)["cj"];
    Rabij["abij"] -= (*Vaibj)["bkci"] * (*Tai)["ak"] * (*Tai)["cj"];
    Rabij["abij"] -= (*Vijka)["jika"] * (*Tai)["bk"];
    Rabij["abij"] += (*Vabij)["acik"] * (*Tai)["cj"] * (*Tai)["bk"];

    // Build Xakic
    Xakic["akic"]  = (*Vabij)["acik"];
    Xakic["akic"] -= (*Vijka)["lkic"] * (*Tai)["al"];
    Xakic["akic"] += (*Vabci)["acdk"] * (*Tai)["di"];
    Xakic["akic"] -= 0.5 * (*Vabij)["dclk"] * (*Tabij)["dail"];
    Xakic["akic"] -= (*Vabij)["dclk"] * (*Tai)["di"] * (*Tai)["al"];
    Xakic["akic"] += (*Vabij)["dclk"] * (*Tabij)["adil"];
    //Xakic["akic"] -= 0.5 * (*Vabij)["cdlk"] * (*Tabij)["adil"]; // Removed in DCSD

    // Build Xakci
    Xakci["akci"]  = (*Vaibj)["akci"];
    Xakci["akci"] -= (*Vijka)["klic"] * (*Tai)["al"];
    Xakci["akci"] += (*Vabci)["adck"] * (*Tai)["di"];
    //Xakci["akci"] -= 0.5 * (*Vabij)["cdlk"] * (*Tabij)["dail"]; // Removed in DCSD
    Xakci["akci"] -= (*Vabij)["cdlk"] * (*Tai)["di"] * (*Tai)["al"];

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
    Xklij["klij"] += (*Vabij)["cdkl"] * (*Tai)["ci"] * (*Tai)["dj"]; 

    // Contract Xklij with T2 Amplitudes
    Rabij["abij"] += Xklij["klij"] * (*Tabij)["abkl"];

    // Contract Xklij with T1 Amplitudes
    Xklij["klij"] += (*Vabij)["cdkl"] * (*Tabij)["cdij"]; //Removed in Dcsd from T2 Amplitudes
    Rabij["abij"] += Xklij["klij"] * (*Tai)["ak"] * (*Tai)["bl"];

    if (Vabcd) {
      // Build Xabcd intermediate
      Tensor<> Xabcd(Vabcd);
      Xabcd["abcd"] -= (*Vabci)["cdak"] * (*Tai)["bk"];
      Xabcd["abcd"] -= (*Vabci)["dcbk"] * (*Tai)["ak"];

      // Contract Xabcd with T2 and T1 Amplitudes
      Rabij["abij"] += Xabcd["abcd"] * (*Tabij)["cdij"];
      Rabij["abij"] += Xabcd["abcd"] * (*Tai)["ci"] * (*Tai)["dj"];
    } else {
      // Slice if Vabcd is not specified

      // Read the sliceRank. If not provided use No
      int64_t sliceRank(getIntegerArgument
			("sliceRank",No));

      // Slice loop starts here
      for (int b(0); b < Nv; b += sliceRank) {
        for (int a(b); a < Nv; a += sliceRank) {
          LOG(0, "DCSD") << "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
          // get the sliced integrals already coupled to the singles
          Tensor<> *Xxycd(sliceCoupledCoulombIntegrals(a, b, sliceRank));
          int lens[] = { Xxycd->lens[0], Xxycd->lens[1], No, No };
          int syms[] = {NS, NS, NS, NS};
          Tensor<> Rxyij(4, lens, syms, *Xxycd->wrld);
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

    LOG(1, "DCSD") << " OK" << std::endl;

    //********************************************************************************
    //***********************  T1 amplitude equations  *******************************
    //********************************************************************************

    LOG(1, "DCSD") << "Solving T1 DCSD Amplitude Equations  ...";

    // Contract Kac and Kki with T1 amplitudes
    Rai["ai"]  = Kac["ac"] * (*Tai)["ci"];
    Rai["ai"] -= Kki["ki"] * (*Tai)["ak"];

    // Build Kck
    Kck["ck"]  = 2.0 * (*Vabij)["cdkl"] * (*Tai)["dl"];
    Kck["ck"] -= (*Vabij)["cdlk"] * (*Tai)["dl"];

    // Contract all the rest terms with T1 and T2 amplitudes
    Rai["ai"] += 2.0 * Kck["ck"] * (*Tabij)["caki"];
    Rai["ai"] -= Kck["ck"] * (*Tabij)["caik"];
    Rai["ai"] += Kck["ck"] * (*Tai)["ci"] * (*Tai)["ak"];
    Rai["ai"] += 2.0 * (*Vabij)["acik"] * (*Tai)["ck"];
    Rai["ai"] -= (*Vaibj)["akci"] * (*Tai)["ck"];
    Rai["ai"] += 2.0 * (*Vabci)["cdak"] * (*Tabij)["cdik"];
    Rai["ai"] -= (*Vabci)["dcak"] * (*Tabij)["cdik"];
    Rai["ai"] += 2.0 * (*Vabci)["cdak"] * (*Tai)["ci"] * (*Tai)["dk"];
    Rai["ai"] -= (*Vabci)["dcak"] * (*Tai)["ci"] * (*Tai)["dk"];
    Rai["ai"] -= 2.0 * (*Vijka)["klic"] * (*Tabij)["ackl"];
    Rai["ai"] += (*Vijka)["lkic"] * (*Tabij)["ackl"];
    Rai["ai"] -= 2.0 * (*Vijka)["klic"] * (*Tai)["ak"] * (*Tai)["cl"];
    Rai["ai"] += (*Vijka)["lkic"] * (*Tai)["ak"] * (*Tai)["cl"];

    singlesAmplitudesFromResiduum(Rai);
    TaiMixer->append(Rai);

    LOG(1, "DCSD") << " OK" << std::endl;
  }
}

