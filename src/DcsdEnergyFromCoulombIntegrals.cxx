#include <DcsdEnergyFromCoulombIntegrals.hpp>
#include <util/Log.hpp>
#include <util/MathFunctions.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <Options.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DcsdEnergyFromCoulombIntegrals);

DcsdEnergyFromCoulombIntegrals::DcsdEnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

DcsdEnergyFromCoulombIntegrals::~DcsdEnergyFromCoulombIntegrals() {
}

/**
 * \brief Calculates DCSD energy from Coulomb integrals Vabcd Vabij Vaibj Vaijb Vijkl Vabci Vijka
 */
void DcsdEnergyFromCoulombIntegrals::run() {
  // Read the Coulomb Integrals Vabij required for the DCSD energy
  Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));

  // Read the Particle/Hole Eigenenergies epsi epsa required for the DCSD energy
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  
  // Compute the no,nv,np
  int no(epsi->lens[0]);
  int nv(epsa->lens[0]);

  // Allocate the DCSD amplitudes Tabij and Tai
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { nv, nv, no, no };
  int vo[]   = { nv, no };
  Tensor<> *Tabij(new Tensor<>(4, vvoo, syms, *Cc4s::world, "Tabij"));
  allocatedTensorArgument("DcsdDoublesAmplitudes", Tabij);
  Tensor<> *Tai(new Tensor<>(2, vo, syms, *Cc4s::world, "Tai"));
  allocatedTensorArgument("DcsdSinglesAmplitudes", Tai);

  // Allocate the DCSD energy e
  Scalar<> energy(*Cc4s::world);
  double e(0), dire, exce;

  LOG(0) <<
    "Solving Distinguishable Cluster Singles and Doubles amplitude equations:" <<
    std::endl;

  // Iteration for determining the DCSD amplitudes Tabij, Tai
  // and the Dcsd energy e
  for (int i(0); i < Cc4s::options->niter; ++i) {
    LOG(0) << "iteration: " << i+1 << std::endl;
    iterate();
    // Singles direct term
    energy[""]  = 4.0 * (*Vabij)["abij"] * (*Tai)["ai"] * (*Tai)["bj"];
    // Doubles direct term
    energy[""] += 2.0 * (*Tabij)["abij"] * (*Vabij)["abij"];
    // Compute direct energy
    dire = energy.get_val();
    // Singles exchange term
    energy[""]  = 2.0 * (*Vabij)["baij"] * (*Tai)["ai"] * (*Tai)["bj"];
    // Doubles exchange term
    energy[""] += (*Tabij)["abji"] * (*Vabij)["abij"];
    // Compute exchange energy
    exce = -1.0 * energy.get_val();
    // Compute total energy
    e = dire + exce;
    LOG(0) << "e=" << e << std::endl;
    LOG(1) << "DCSDdir=" << dire << std::endl;
    LOG(1) << "DCSDexc=" << exce << std::endl;
  }

  LOG(1) << "DCSD correlation energy = " << e << std::endl;

  setRealArgument("DcsdEnergy", e);
}

//////////////////////////////////////////////////////////////////////
// Hirata iteration routine for the CCSD amplitudes Tabij and Tai from
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
// modified to give DCSD amplitudes according to
// D. Kats, et. al., J. Chem. Phys. 142, 064111 (2015)
//////////////////////////////////////////////////////////////////////
void DcsdEnergyFromCoulombIntegrals::iterate() {
  {
    // Read the DCSD amplitudes Tai and Tabij
    Tensor<> *Tai(getTensorArgument("DcsdSinglesAmplitudes"));
    Tensor<> *Tabij(getTensorArgument("DcsdDoublesAmplitudes"));

    // Read the Coulomb Integrals Vabcd Vabij Vaibj Vijkl Vabci Vijka
    Tensor<> *Vabcd(getTensorArgument("PPPPCoulombIntegrals"));
    Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
    Tensor<> *Vaibj(getTensorArgument("PHPHCoulombIntegrals"));
    Tensor<> *Vijkl(getTensorArgument("HHHHCoulombIntegrals"));
    Tensor<> *Vijka(getTensorArgument("HHHPCoulombIntegrals"));
    Tensor<> *Vabci(getTensorArgument("PPPHCoulombIntegrals"));

    // Read the Particle/Hole Eigenenergies epsi epsa
    Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
    Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  
    // Compute the no,nv,np
    int no(epsi->lens[0]);
    int nv(epsa->lens[0]);

    int syms[] = { NS, NS, NS, NS };
    int voov[] = { nv, no, no, nv };
    int vo[] = { nv, no };
    int vv[] = { nv, nv };
    int oo[] = { no, no };

    // Allocate Tensors for T1 amplitudes
    Tensor<> Rai(false, *Tai);
    Tensor<> Dai(false, *Tai);

    // Allocate Tensors for T2 amplitudes
    Tensor<> Rabij(false, *Vabij);
    Tensor<> Dabij(false, *Vabij);

    // Define intermediates
    Tensor<> Lac(2, vv, syms, *Cc4s::world, "Lac");
    Tensor<> Kac(2, vv, syms, *Cc4s::world, "Kac");
    Tensor<> Lki(2, oo, syms, *Cc4s::world, "Lki");
    Tensor<> Kki(2, oo, syms, *Cc4s::world, "Kki");
    Tensor<> Kck(2, vo, syms, *Cc4s::world, "Kck");

    Tensor<> Xklij(false, *Vijkl);
    Tensor<> Xakci(false, *Vaibj);
    Tensor<> Xabcd(false, *Vabcd);
    Tensor<> Xakic(4, voov, syms, *Cc4s::world, "Xakic");

    //********************************************************************************
    //***********************  T2 amplitude equations  *******************************
    //********************************************************************************

    LOG(1) << "Solving T2 DCSD Amplitude Equations  ...";

    // Build Kac
    Kac["ac"]  = -2.0 * (*Vabij)["cdkl"] * (*Tabij)["adkl"];
    Kac["ac"] += (*Vabij)["dckl"] * (*Tabij)["adkl"];
    Kac["ac"] -= 2.0 * (*Vabij)["cdkl"] * (*Tai)["ak"] * (*Tai)["dl"];
    Kac["ac"] += (*Vabij)["dckl"] * (*Tai)["ak"] * (*Tai)["dl"];

    // Build Lac
    Lac["ac"]  = 0.5 * Kac["ac"]; // Multiplied by 0.5 in DCSD
    Lac["ac"] += 2.0 * (*Vabci)["cdak"] * (*Tai)["dk"];
    Lac["ac"] -= (*Vabci)["dcak"] * (*Tai)["dk"];

    // Build Kki
    Kki["ki"]  = 2.0 * (*Vabij)["cdkl"] * (*Tabij)["cdil"];
    Kki["ki"] -= (*Vabij)["dckl"] * (*Tabij)["cdil"];
    Kki["ki"] += 2.0 * (*Vabij)["cdkl"] * (*Tai)["ci"] * (*Tai)["dl"];
    Kki["ki"] -= (*Vabij)["dckl"] * (*Tai)["ci"] * (*Tai)["dl"];

    // Build Lki
    Lki["ki"]  = 0.5 * Kki["ki"]; // Multiplied by 0.5 in DCSD
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
    Xakic["akic"] += (*Vabci)["adck"] * (*Tai)["di"];
    Xakic["akic"] -= 0.5 * (*Vabij)["dclk"] * (*Tabij)["dail"];
    Xakic["akic"] -= (*Vabij)["dclk"] * (*Tai)["di"] * (*Tai)["al"];
    Xakic["akic"] += (*Vabij)["dclk"] * (*Tabij)["adil"];
    //Xakic["akic"] -= 0.5 * (*Vabij)["cdlk"] * (*Tabij)["adil"]; // Removed in DCD

    // Build Xakci
    Xakci["akci"]  = (*Vaibj)["akci"];
    Xakci["akci"] -= (*Vijka)["klic"] * (*Tai)["al"];
    Xakci["akci"] += (*Vabci)["adck"] * (*Tai)["di"];
    //Xakci["akci"] -= 0.5 * (*Vabij)["cdlk"] * (*Tabij)["dail"]; // Removed in DCD
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
    //Xklij["klij"] += (*Vabij)["cdkl"] * (*Tabij)["cdij"]; //Removed in Dcd
    Xklij["klij"] += (*Vabij)["cdkl"] * (*Tai)["ci"] * (*Tai)["dj"]; 

    // Contract Xklij with T2 Amplitudes
    Rabij["abij"] += Xklij["klij"] * (*Tabij)["abkl"];

    // Contract Xklij with T1 Amplitudes
    Rabij["abij"] += Xklij["klij"] * (*Tai)["ak"] * (*Tai)["bl"];

    // Build Xabcd intermediate
    Xabcd["abcd"]  = (*Vabcd)["abcd"];
    Xabcd["abcd"] -= (*Vabci)["cdak"] * (*Tai)["bk"];
    Xabcd["abcd"] -= (*Vabci)["dcbk"] * (*Tai)["ak"];
    
    // Contract Xabcd with T2 and T1 Amplitudes
    Rabij["abij"] += Xabcd["abcd"] * (*Tabij)["cdij"];
    Rabij["abij"] += Xabcd["abcd"] * (*Tai)["ci"] * (*Tai)["dj"];

    // Build Dabij
    Dabij["abij"]  = (*epsi)["i"];
    Dabij["abij"] += (*epsi)["j"];
    Dabij["abij"] -= (*epsa)["a"];
    Dabij["abij"] -= (*epsa)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
    Dabij["abij"] = Dabij["abij"];

    // Divide Rabij/Dabij to get Tabij
    Bivar_Function<> fDivide(&divide<double>);
    Tabij->contract(1.0, Rabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);

    LOG(1) << " OK" << std::endl;

    //********************************************************************************
    //***********************  T1 amplitude equations  *******************************
    //********************************************************************************

    LOG(1) << "Solving T1 DCSD Amplitude Equations  ...";

    // Contract Kac and Kki with T1 amplitudes
    Rai["ai"]  = Kac["ac"] * (*Tai)["ci"];
    Rai["ai"] -= Kki["ki"] * (*Tai)["ak"];

    //Build Kck
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

    // Build Dai
    Dai["ai"]  = (*epsi)["i"];
    Dai["ai"] -= (*epsa)["a"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
    Dai["ai"] = Dai["ai"];

    // Divide Rai/Dai to get Tai
    Tai->contract(1.0, Rai,"ai", Dai,"ai", 0.0,"ai", fDivide);

    LOG(1) << " OK" << std::endl;
  }
}
