#include <algorithms/DcdEnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <Options.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DcdEnergyFromCoulombIntegrals);

DcdEnergyFromCoulombIntegrals::DcdEnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

DcdEnergyFromCoulombIntegrals::~DcdEnergyFromCoulombIntegrals() {
}

/**
 * \brief Calculates DCD energy from Coulomb integrals Vabcd Vabij Vaibj Vijkl
 */

void DcdEnergyFromCoulombIntegrals::run() {
  // Read the Coulomb Integrals Vabij required for the DCD energy
  Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));

  // Read the Particle/Hole Eigenenergies epsi epsa required for the DCD energy
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  
  // Compute the no,nv,np
  int no(epsi->lens[0]);
  int nv(epsa->lens[0]);

  // Allocate the DCD amplitudes Tabij
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { nv, nv, no, no };
  Tensor<> *Tabij(new Tensor<>(4, vvoo, syms, *Cc4s::world, "Tabij"));
  allocatedTensorArgument("DcdDoublesAmplitudes", Tabij);

  // Allocate the DCD energy e
  Scalar<> energy(*Cc4s::world);
  double e(0), dire, exce;

  LOG(0) <<
    "Solving Distinguishable Cluster Doubles Amplitude Equations:" <<
    std::endl;

  // Iteration for determining the DCD amplitudes Tabij
  // and the Dcd energy e
  for (int i(0); i < Cc4s::options->niter; ++i) {
    LOG(0) << "iteration: " << i+1 << std::endl;
    iterateHirata();
    energy[""] = 2.0 * (*Tabij)["abij"] * (*Vabij)["abij"];
    dire = energy.get_val();
    energy[""] = (*Tabij)["abji"] * (*Vabij)["abij"];
    exce = -1.0 * energy.get_val();
    e = dire + exce;
    LOG(0) << "e=" << e << std::endl;
    LOG(1) << "DCDdir=" << dire << std::endl;
    LOG(1) << "DCDexc=" << exce << std::endl;
    // Print the MP2 energy in 1st iteration
    if (i == 0) {
      LOG(1) << "MP2 correlation energy = " << e << std::endl;      
    }
  }

  LOG(1) << "DCD correlation energy = " << e << std::endl;

  setRealArgument("DcdEnergy", e);
}

//////////////////////////////////////////////////////////////////////
// Hiarata iteration routine for the DCD amplitudes Tabij (Table. 1)
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
// Modified according to D. Kats, J. Chem. Phys. 139, 021102 (2013)
//////////////////////////////////////////////////////////////////////
void DcdEnergyFromCoulombIntegrals::iterateHirata() {
  {
    // Read the DCD amplitudes Tabij
    Tensor<> *Tabij(getTensorArgument("DcdDoublesAmplitudes"));

    // Read the Coulomb Integrals Vabcd Vabij Vaibj Vijkl
    Tensor<> *Vabcd(getTensorArgument("PPPPCoulombIntegrals"));
    Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
    Tensor<> *Vaibj(getTensorArgument("PHPHCoulombIntegrals"));
    Tensor<> *Vijkl(getTensorArgument("HHHHCoulombIntegrals"));

    // Read the Particle/Hole Eigenenergies epsi epsa
    Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
    Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  
    // Compute the no,nv,np
    int no(epsi->lens[0]);
    int nv(epsa->lens[0]);

    int syms[] = { NS, NS, NS, NS };
    int voov[] = { nv, no, no, nv };
    int vv[] = { nv, nv };
    int oo[] = { no, no };

    // Allocate Tensors for T2 amplitudes
    Tensor<> Rabij(false, *Vabij);
    Tensor<> Dabij(false, *Vabij);

    // Define intermediates
    Tensor<> Kac(2, vv, syms, *Cc4s::world, "Kac");
    Tensor<> Kki(2, oo, syms, *Cc4s::world, "Kki");

    Tensor<> Xklij(false, *Vijkl);
    Tensor<> Xakci(false, *Vaibj);
    Tensor<> Xakic(4, voov, syms, *Cc4s::world, "Xakic");

    // Build Kac
    Kac["ac"]  = -2.0 * (*Vabij)["cdkl"] * (*Tabij)["adkl"];
    Kac["ac"] += (*Vabij)["dckl"] * (*Tabij)["adkl"];

    // Build Kki
    Kki["ki"]  = 2.0 * (*Vabij)["cdkl"] * (*Tabij)["cdil"];
    Kki["ki"] -= (*Vabij)["dckl"] * (*Tabij)["cdil"];
    
    // Contract Kac with T2 Amplitudes
    Rabij["abij"]  = 0.5 * Kac["ac"] * (*Tabij)["cbij"]; // Multiplied by 0.5 in DCD

    // Contract Kki with T2 Amplitudes
    Rabij["abij"] -= 0.5 * Kki["ki"] * (*Tabij)["abkj"]; // Multiplied by 0.5 in DCD

    // Build Xakic
    Xakic["akic"]  = (*Vabij)["acik"];
    Xakic["akic"] -= 0.5 * (*Vabij)["dclk"] * (*Tabij)["dail"];
    Xakic["akic"] += (*Vabij)["dclk"] * (*Tabij)["adil"];
    //Xakic["akic"] -= 0.5 * (*Vabij)["cdlk"] * (*Tabij)["adil"]; // Removed in DCD

    // Build Xakci
    Xakci["akci"]  = (*Vaibj)["akci"];
    //Xakci["akci"] -= 0.5 * (*Vabij)["cdlk"] * (*Tabij)["dail"]; // Removed in DCD

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
    // Xklij["klij"] += (*Vabij)["cdkl"] * (*Tabij)["cdij"]; //Removed in Dcd

    // Contract Xklij with T2 Amplitudes
    Rabij["abij"] += Xklij["klij"] * (*Tabij)["abkl"];

    // Contract Vabcd with T2 Amplitudes
    Rabij["abij"] += (*Vabcd)["abcd"] * (*Tabij)["cdij"];

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
  }
}

//////////////////////////////////////////////////////////////////////
// Bartlett iteration routine for the DCD amplitudes Tabij 
// Rev. Mod. Phys. 79, 291  Page 305, Figure 8. -> CCD
// J. Chem. Phys. 139, 021102 (2013) -> DCD
//////////////////////////////////////////////////////////////////////
void DcdEnergyFromCoulombIntegrals::iterateBartlett() {
  {
    // Read the DCD amplitudes Tabij
    Tensor<> *Tabij(getTensorArgument("DcdDoublesAmplitudes"));

    // Read the Coulomb Integrals Vabcd Vabij Vaibj Vijkl
    Tensor<> *Vabcd(getTensorArgument("PPPPCoulombIntegrals"));
    Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
    Tensor<> *Vaibj(getTensorArgument("PHPHCoulombIntegrals"));
    Tensor<> *Vijkl(getTensorArgument("HHHHCoulombIntegrals"));

    // Read the Particle/Hole Eigenenergies epsi epsa
    Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
    Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  
    // Compute the no,nv,np
    int no(epsi->lens[0]);
    int nv(epsa->lens[0]);

    int syms[] = { NS, NS, NS, NS };
    int vv[] = { nv, nv };
    int oo[] = { no, no };

    // Allocate Tensors for T2 amplitudes
    Tensor<> Rabij(false, *Vabij);
    Tensor<> Dabij(false, *Vabij);

    // Define intermediates
    Tensor<> Cabij(false, *Vabij);
    Tensor<> Caibj(false, *Vaibj);
    Tensor<> Cabcd(false, *Vabcd);
    Tensor<> Cij(2, oo, syms, *Cc4s::world, "Cij");
    Tensor<> Cab(2, vv, syms, *Cc4s::world, "Cab");

    //////////////////////////////////////////////////////////////////////
    // Create linear terms with T2 Amplitudes that need permutation
    //////////////////////////////////////////////////////////////////////
    // Contract Vabcd with T2 Amplitudes (3rd term first line)
    Rabij["abij"]  = 0.5*(*Vabcd)["abcd"] * (*Tabij)["cdij"];

    // Contract Vijkl with T2 Amplitudes (4th term first line)
    Rabij["abij"] += 0.5*(*Vijkl)["klij"] * (*Tabij)["abkl"];

    // Contract Vabij with T2 Amplitudes (1st term second line)
    Rabij["abij"] += 2.0*(*Vabij)["cbkj"] * (*Tabij)["acik"];

    // Contract Vabij with T2 Amplitudes (2nd term second line)
    Rabij["abij"] -= 1.0*(*Vabij)["cbkj"] * (*Tabij)["caik"];

    // Contract Vaibj with T2 Amplitudes (3rd term second line)
    Rabij["abij"] -= 1.0*(*Vaibj)["cibk"] * (*Tabij)["ackj"];

    // Contract Vaibj with T2 Amplitudes (4th term second line)
    Rabij["abij"] -= 1.0*(*Vaibj)["cjbk"] * (*Tabij)["acik"];


    //////////////////////////////////////////////////////////////////////
    // Create quadratic terms with T2 Amplitudes that need permutation
    //////////////////////////////////////////////////////////////////////
    // Contract Vabij with T2*T2 Amplitudes
    //////////////////////////////////////////////////////////////////////
    // 1st term third line
    Cabij["dblj"]  = 1.0*(*Vabij)["dclk"] * (*Tabij)["cbkj"];
    Rabij["abij"] += 2.0*   Cabij["cbkj"] * (*Tabij)["acik"];

    // 2nd term third line
    Cabij["dblj"]  = 1.0*(*Vabij)["dclk"] * (*Tabij)["cbjk"];
    Rabij["abij"] -= 2.0*   Cabij["cbkj"] * (*Tabij)["acik"];

    // 3rd term third line
    Cabij["dblj"]  = 1.0*(*Vabij)["dclk"] * (*Tabij)["cbjk"];
    Rabij["abij"] += 0.5*   Cabij["cbkj"] * (*Tabij)["caik"];

    // 1st term fourth line (Removed in Dcd)
    // Cabij["dblj"]  = 1.0*(*Vabij)["cdlk"] * (*Tabij)["cbkj"];
    // Rabij["abij"] -= 1.0*   Cabij["dblj"] * (*Tabij)["adil"];

    // 2nd term fourth line (Removed in Dcd)
    // Cabij["dblj"]  = 1.0*(*Vabij)["cdlk"] * (*Tabij)["cbkj"];
    // Rabij["abij"] += 1.0*   Cabij["dblj"] * (*Tabij)["adli"];

    // 3rd term fourth line (Removed in Dcd)
    // Cabij["dbli"]  = 1.0*(*Vabij)["cdlk"] * (*Tabij)["cbik"];
    // Rabij["abij"] += 0.5*   Cabij["dbli"] * (*Tabij)["adlj"];

    // 1st term fifth line (Removed in Dcd)
    // Cabcd["cdab"]  = 1.0*(*Vabij)["cdlk"] * (*Tabij)["ablk"];
    // Rabij["abij"] += 0.5*   Cabcd["cdab"] * (*Tabij)["cdij"];

    // 2nd term fifth line (Multiplied by *0.5 in Dcd)
    Cij["li"]      = 1.0*(*Vabij)["cdkl"] * (*Tabij)["cdki"];
    Rabij["abij"] -= 1.0*       Cij["li"] * (*Tabij)["ablj"];

    // 3rd term fifth line (Multiplied by *0.5 in Dcd)
    Cij["li"]      = 1.0*(*Vabij)["cdkl"] * (*Tabij)["cdik"];
    Rabij["abij"] += 0.5*       Cij["li"] * (*Tabij)["ablj"];

    // 1st term sixth line (Multiplied by *0.5 in Dcd)
    Cab["da"]      = 1.0*(*Vabij)["cdkl"] * (*Tabij)["cakl"];
    Rabij["abij"] -= 1.0*       Cab["da"] * (*Tabij)["dbij"];

    // 2nd term sixth line (Multiplied by *0.5 in Dcd)
    Cab["da"]      = 1.0*(*Vabij)["cdkl"] * (*Tabij)["ackl"];
    Rabij["abij"] += 0.5*       Cab["da"] * (*Tabij)["dbij"];

    //////////////////////////////////////////////////////////////////////
    // Symmetrize Rabij by applying permutation operator
    // To save memory we use Caibj as intermediate for the permutation operator 
    //////////////////////////////////////////////////////////////////////
    Caibj["aibj"]  = Rabij["abij"];
    Rabij["abij"] += Caibj["bjai"];

    //////////////////////////////////////////////////////////////////////
    // Add the Vabij (the only term that does not need permutation)
    //////////////////////////////////////////////////////////////////////
    Rabij["abij"] += (*Vabij)["abij"];

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
  }
}
