#include <CcdEnergyFromCoulombIntegrals.hpp>
#include <util/Log.hpp>
#include <util/MathFunctions.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <Options.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CcdEnergyFromCoulombIntegrals);

CcdEnergyFromCoulombIntegrals::CcdEnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

CcdEnergyFromCoulombIntegrals::~CcdEnergyFromCoulombIntegrals() {
}

/**
 * \brief Calculates CCD energy from Coulomb integrals Vabcd Vabij Vaibj Vijkl
 */
void CcdEnergyFromCoulombIntegrals::run() {
  // Read the Coulomb Integrals Vabij required for the CCD energy
  Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));

  // Read the Particle/Hole Eigenenergies epsi epsa required for the CCD energy
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  
  // Compute the no,nv,np
  int no(epsi->lens[0]);
  int nv(epsa->lens[0]);

  // Allocate the CCD amplitudes tabij
  int syms[] = { NS, NS, NS, NS };
  int vvoo[] = { nv, nv, no, no };
  Tensor<> *Tabij(new Tensor<>(4, vvoo, syms, *Cc4s::world, "Tabij"));
  allocatedTensorArgument("CcdDoublesAmplitudes", Tabij);

  // Allocate the CCD energy e
  Scalar<> energy(*Cc4s::world);
  double e(0), dire, exce;

  LOG(0) <<
    "Solving Coupled Cluster Doubles Amplitude Equations:" <<
    std::endl;

  // Iteration for determining the CCD amplitudes Tabij
  // and the Ccd energy e
  for (int i(0); i < Cc4s::options->niter; ++i) {
    LOG(0) << "iteration: " << i+1 << std::endl;
    iterateHirata();
    energy[""] = 2.0 * (*Tabij)["abij"] * (*Vabij)["abij"];
    dire = energy.get_val();
    energy[""] = (*Tabij)["abji"] * (*Vabij)["abij"];
    exce = -1.0 * energy.get_val();
    e = dire + exce;
    LOG(0) << "e=" << e << std::endl;
    LOG(1) << "CCDdir=" << dire << std::endl;
    LOG(1) << "CCDexc=" << exce << std::endl;
  }

  LOG(0) << "CCD correlation energy = " << e << std::endl;

  setRealArgument("CcdEnergy", e);
}

//////////////////////////////////////////////////////////////////////
// Hiarata iteration routine for the CCD amplitudes Tabij (Table. 1)
// So Hirata, et. al. Chem Phys Letters, 345, 475 (2001)
// http://dx.doi.org/10.1016/S0009-2614(01)00897-1
//////////////////////////////////////////////////////////////////////
void CcdEnergyFromCoulombIntegrals::iterateHirata() {
  {
    // Read the CCD amplitudes Tabij
    Tensor<> *Tabij(getTensorArgument("CcdDoublesAmplitudes"));

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
    int vo[] = { nv, no };
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
    Rabij["abij"] = Kac["ac"] * (*Tabij)["cbij"];

    // Contract Kki with T2 Amplitudes
    Rabij["abij"] -= Kki["ki"] * (*Tabij)["abkj"];

    // Build Xakic
    Xakic["akic"]  = (*Vabij)["acik"];
    Xakic["akic"] -= 0.5 * (*Vabij)["dclk"] * (*Tabij)["dail"];
    Xakic["akic"] += (*Vabij)["dclk"] * (*Tabij)["adil"];
    Xakic["akic"] -= 0.5 * (*Vabij)["cdlk"] * (*Tabij)["adil"];

    // Build Xakci
    Xakci["akci"]  = (*Vaibj)["akci"];
    Xakci["akci"] -= 0.5 * (*Vabij)["cdlk"] * (*Tabij)["dail"];

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
    Xklij["klij"] += (*Vabij)["cdkl"] * (*Tabij)["cdij"];

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
// Bartlett iteration routine for the CCD amplitudes Tabij 
// (Rev. Mod. Phys. 79, 291  Page 305, Figure 8.)
//////////////////////////////////////////////////////////////////////
void CcdEnergyFromCoulombIntegrals::iterateBartlett() {
  {
    // Read the CCD amplitudes Tabij
    Tensor<> *Tabij(getTensorArgument("CcdDoublesAmplitudes"));

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
    // 1st term third line
    Cabij["dblj"]  = 1.0*(*Vabij)["dclk"] * (*Tabij)["cbkj"];
    Rabij["abij"] += 2.0*   Cabij["cbkj"] * (*Tabij)["acik"];

    // 2nd term third line
    Cabij["dblj"]  = 1.0*(*Vabij)["dclk"] * (*Tabij)["cbjk"];
    Rabij["abij"] -= 2.0*   Cabij["cbkj"] * (*Tabij)["acik"];

    // 3rd term third line
    Cabij["dblj"]  = 1.0*(*Vabij)["dclk"] * (*Tabij)["cbjk"];
    Rabij["abij"] += 0.5*   Cabij["cbkj"] * (*Tabij)["caik"];

    // 1st term fourth line
    Cabij["dblj"]  = 1.0*(*Vabij)["cdlk"] * (*Tabij)["cbkj"];
    Rabij["abij"] -= 1.0*   Cabij["dblj"] * (*Tabij)["adil"];

    // 2nd term fourth line
    Cabij["dblj"]  = 1.0*(*Vabij)["cdlk"] * (*Tabij)["cbkj"];
    Rabij["abij"] += 1.0*   Cabij["dblj"] * (*Tabij)["adli"];

    // 3rd term fourth line
    Cabij["dbli"]  = 1.0*(*Vabij)["cdlk"] * (*Tabij)["cbik"];
    Rabij["abij"] += 0.5*   Cabij["dbli"] * (*Tabij)["adlj"];

    // 1st term fifth line
    Cabcd["cdab"]  = 1.0*(*Vabij)["cdlk"] * (*Tabij)["ablk"];
    Rabij["abij"] += 0.5*   Cabcd["cdab"] * (*Tabij)["cdij"];

    // 2nd term fifth line
    Cij["li"]      = 1.0*(*Vabij)["cdkl"] * (*Tabij)["cdki"];
    Rabij["abij"] -= 2.0*       Cij["li"] * (*Tabij)["ablj"];

    // 3rd term fifth line
    Cij["li"]      = 1.0*(*Vabij)["cdkl"] * (*Tabij)["cdik"];
    Rabij["abij"] += 1.0*       Cij["li"] * (*Tabij)["ablj"];

    // 1st term sixth line
    Cab["da"]      = 1.0*(*Vabij)["cdkl"] * (*Tabij)["cakl"];
    Rabij["abij"] -= 2.0*       Cab["da"] * (*Tabij)["dbij"];

    // 2nd term sixth line
    Cab["da"]      = 1.0*(*Vabij)["cdkl"] * (*Tabij)["ackl"];
    Rabij["abij"] += 1.0*       Cab["da"] * (*Tabij)["dbij"];

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
