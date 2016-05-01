#include <algorithms/DcdEnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <util/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DcdEnergyFromCoulombIntegrals);

DcdEnergyFromCoulombIntegrals::DcdEnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): ClusterDoublesAlgorithm(argumentList) {
}

DcdEnergyFromCoulombIntegrals::~DcdEnergyFromCoulombIntegrals() {
}

//////////////////////////////////////////////////////////////////////
// Hiarata iteration routine for the DCD amplitudes Tabij (Table. 1)
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
// Modified according to D. Kats, J. Chem. Phys. 139, 021102 (2013)
//////////////////////////////////////////////////////////////////////
void DcdEnergyFromCoulombIntegrals::iterate(int i) {
  {
    // Read the DCD amplitudes Tabij
    Tensor<> *Tabij(&TabijMixer->getNext());

    // Read the Coulomb Integrals Vabij
    Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));

    // Allocate Tensor for T2 amplitudes
    Tensor<> Rabij(false, *Vabij);

    if (i == 0) {
      // For first iteration compute only the MP2 amplitudes since Tabij = 0

      // Rabij contracted with Vabij is the only non-zero term
      Rabij["abij"] += (*Vabij)["abij"];
      
    } else {
      // For the rest iterations compute the DCD amplitudes

      // Read the Coulomb Integrals Vabcd Vaibj Vijkl
      // the PPPPCoulombIntegrals may not be given then slicing is required
      Tensor<> *Vabcd(isArgumentGiven("PPPPCoulombIntegrals") ?
		      getTensorArgument("PPPPCoulombIntegrals") : nullptr);
      Tensor<> *Vaibj(getTensorArgument("PHPHCoulombIntegrals"));
      Tensor<> *Vijkl(getTensorArgument("HHHHCoulombIntegrals"));

      // Read the Particle/Hole Eigenenergies epsi epsa
      Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
      Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
      
      // Compute the No,Nv
      int No(epsi->lens[0]);
      int Nv(epsa->lens[0]);

      std::string abbreviation(getAbbreviation());
      std::transform(abbreviation.begin(), abbreviation.end(), 
		     abbreviation.begin(), ::toupper);

      {
	// Define intermediates
	int syms[] = { NS, NS, NS, NS };
	int voov[] = { Nv, No, No, Nv };
	int vv[] = { Nv, Nv };
	int oo[] = { No, No };

	Tensor<> Kac(2, vv, syms, *epsi->wrld, "Kac");
	Tensor<> Kki(2, oo, syms, *epsi->wrld, "Kki");

	Tensor<> Xklij(false, *Vijkl);
	Tensor<> Xakci(false, *Vaibj);
	Tensor<> Xakic(4, voov, syms, *epsi->wrld, "Xakic");

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
	// Xklij["klij"] += (*Vabij)["cdkl"] * (*Tabij)["cdij"]; //Removed in DCD

	// Contract Xklij with T2 Amplitudes
	Rabij["abij"] += Xklij["klij"] * (*Tabij)["abkl"];
      }
      
      // Contract Vabcd with T2 Amplitudes
      if (Vabcd) {
	Rabij["abij"] += (*Vabcd)["abcd"] * (*Tabij)["cdij"];
      } else {
	// Slice if Vabcd is not specified

	// Read the sliceRank. If not provided use No
	int sliceRank(getIntegerArgument
		      ("sliceRank",No));

	// Slice loop starts here
	for (int b(0); b < Nv; b += sliceRank) {
	  for (int a(b); a < Nv; a += sliceRank) {
	    LOG(1, abbreviation) << "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
	    Tensor<> *Vxycd(sliceCoulombIntegrals(a, b, sliceRank));
	    int lens[] = { Vxycd->lens[0], Vxycd->lens[1], No, No };
	    int syms[] = {NS, NS, NS, NS};
	    Tensor<> Rxyij(4, lens, syms, *Vxycd->wrld);
	    Rxyij["xyij"] = (*Vxycd)["xycd"] * (*Tabij)["cdij"];
	    sliceIntoResiduum(Rxyij, a, b, Rabij);
	    // the integrals of this slice are not needed anymore
	    delete Vxycd;
	  }
	}
      }
    }

    // calculate the amplitdues from the residuum
    doublesAmplitudesFromResiduum(Rabij);
    // and append them to the mixer
    TabijMixer->append(Rabij);
  }
}

void DcdEnergyFromCoulombIntegrals::dryIterate() {
  {
    // TODO: the Mixer should provide a DryTensor in the future
    // Read the DCD amplitudes Tabij
    //DryTensor<> *Tabij(
    getTensorArgument<double, DryTensor<double>>("DcdDoublesAmplitudes");
    //);

    // Read the Coulomb Integrals Vabcd Vabij Vaibj Vijkl
    // the PPPPCoulombIntegrals may not be given then slicing is required
    DryTensor<> *Vabcd(
      isArgumentGiven("PPPPCoulombIntegrals") ?
        getTensorArgument<double, DryTensor<double>>("PPPPCoulombIntegrals") :
        nullptr
    );
    DryTensor<> *Vabij(getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals"));
    DryTensor<> *Vaibj(getTensorArgument<double, DryTensor<double>>("PHPHCoulombIntegrals"));
    DryTensor<> *Vijkl(getTensorArgument<double, DryTensor<double>>("HHHHCoulombIntegrals"));

    // Read the Particle/Hole Eigenenergies epsi epsa
    DryTensor<> *epsi(getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies"));
    DryTensor<> *epsa(getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies"));
  
    // Compute the no,nv,np
    int No(epsi->lens[0]);
    int Nv(epsa->lens[0]);

    int syms[] = { NS, NS, NS, NS };
    int voov[] = { Nv, No, No, Nv };
    int vv[] = { Nv, Nv };
    int oo[] = { No, No };

    // Allocate Tensors for T2 amplitudes
    DryTensor<> Rabij(*Vabij);

    // Define intermediates
    DryTensor<> Kac(2, vv, syms);
    DryTensor<> Kki(2, oo, syms);

    DryTensor<> Xklij(*Vijkl);
    DryTensor<> Xakci(*Vaibj);
    DryTensor<> Xakic(4, voov, syms);

    // Read the sliceRank. If not provided use No
    int sliceRank(getIntegerArgument
		  ("sliceRank",No));

    if (!Vabcd) {
      // slice if Vabcd is not specified
      int lens[] = { sliceRank, sliceRank, Nv, Nv };
      int syms[] = {NS, NS, NS, NS};
      // TODO: implement drySliceCoulombIntegrals
      DryTensor<> Vxycd(4, lens, syms);
      DryTensor<> Rxyij(*Vijkl);
    }
    // TODO: implment dryDoublesAmplitudesFromResiduum
    // at the moment, assume usage of Dabij
    DryTensor<> Dabij(*Vabij);
  }
}

//////////////////////////////////////////////////////////////////////
// Bartlett iteration routine for the DCD amplitudes Tabij 
// Rev. Mod. Phys. 79, 291  Page 305, Figure 8. -> CCD
// J. Chem. Phys. 139, 021102 (2013) -> DCD
//////////////////////////////////////////////////////////////////////
void DcdEnergyFromCoulombIntegrals::iterateBartlett(int i) {
  {
    // Read the DCD amplitudes Tabij
    Tensor<> *Tabij(&TabijMixer->getNext());

    // Read the Coulomb Integrals Vabcd Vabij Vaibj Vijkl
    // the PPPPCoulombIntegrals may not be given then slicing is required
    Tensor<> *Vabcd(isArgumentGiven("PPPPCoulombIntegrals") ?
		    getTensorArgument("PPPPCoulombIntegrals") : nullptr
		    );
    Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
    Tensor<> *Vaibj(getTensorArgument("PHPHCoulombIntegrals"));
    Tensor<> *Vijkl(getTensorArgument("HHHHCoulombIntegrals"));

    // Read the Particle/Hole Eigenenergies epsi epsa
    Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
    Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  
    // Compute the No,Nv
    int No(epsi->lens[0]);
    int Nv(epsa->lens[0]);

    int syms[] = { NS, NS, NS, NS };
    int vv[] = { Nv, Nv };
    int oo[] = { No, No };

    std::string abbreviation(getAbbreviation());
  
    // Allocate Tensors for T2 amplitudes
    Tensor<> Rabij(false, *Vabij);

    // Define intermediates
    Tensor<> Cabij(false, *Vabij);
    Tensor<> Caibj(false, *Vaibj);
    // Tensor<> Cabcd(false, *Vabcd);  Not needed in DCD
    Tensor<> Cij(2, oo, syms, *epsi->wrld, "Cij");
    Tensor<> Cab(2, vv, syms, *epsi->wrld, "Cab");

    //////////////////////////////////////////////////////////////////////
    // Create linear terms with T2 Amplitudes that need permutation
    //////////////////////////////////////////////////////////////////////

    // Contract Vijkl with T2 Amplitudes (4th term first line)
    Rabij["abij"] = 0.5*(*Vijkl)["klij"] * (*Tabij)["abkl"];

    // Contract Vabcd with T2 Amplitudes (3rd term first line)
    if (Vabcd) {
      Rabij["abij"] += 0.5*(*Vabcd)["abcd"] * (*Tabij)["cdij"];
    } else {
      // Slice if Vabcd is not specified

      // Read the sliceRank. If not provided use No
      int64_t sliceRank(getIntegerArgument
			("sliceRank",No));

      // Slice loop starts here
      for (int b(0); b < Nv; b += sliceRank) {
        for (int a(b); a < Nv; a += sliceRank) {
          LOG(1, abbreviation) << "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
          Tensor<> *Vxycd(sliceCoulombIntegrals(a, b, sliceRank));
          int lens[] = { Vxycd->lens[0], Vxycd->lens[1], No, No };
          int syms[] = {NS, NS, NS, NS};
          Tensor<> Rxyij(4, lens, syms, *Vxycd->wrld);
          Rxyij["xyij"] = 0.5*(*Vxycd)["xycd"] * (*Tabij)["cdij"];
          sliceIntoResiduum(Rxyij, a, b, Rabij);
          // the integrals of this slice are not needed anymore
          delete Vxycd;
        }
      }
    }

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

    // 1st term fourth line (Removed in DCD)
    // Cabij["dblj"]  = 1.0*(*Vabij)["cdlk"] * (*Tabij)["cbkj"];
    // Rabij["abij"] -= 1.0*   Cabij["dblj"] * (*Tabij)["adil"];

    // 2nd term fourth line (Removed in DCD)
    // Cabij["dblj"]  = 1.0*(*Vabij)["cdlk"] * (*Tabij)["cbkj"];
    // Rabij["abij"] += 1.0*   Cabij["dblj"] * (*Tabij)["adli"];

    // 3rd term fourth line (Removed in DCD)
    // Cabij["dbli"]  = 1.0*(*Vabij)["cdlk"] * (*Tabij)["cbik"];
    // Rabij["abij"] += 0.5*   Cabij["dbli"] * (*Tabij)["adlj"];

    // 1st term fifth line (Removed in DCD)
    // Cabcd["cdab"]  = 1.0*(*Vabij)["cdlk"] * (*Tabij)["ablk"];
    // Rabij["abij"] += 0.5*   Cabcd["cdab"] * (*Tabij)["cdij"];

    // 2nd term fifth line (Multiplied by *0.5 in DCD)
    Cij["li"]      = 1.0*(*Vabij)["cdkl"] * (*Tabij)["cdki"];
    Rabij["abij"] -= 1.0*       Cij["li"] * (*Tabij)["ablj"];

    // 3rd term fifth line (Multiplied by *0.5 in DCD)
    Cij["li"]      = 1.0*(*Vabij)["cdkl"] * (*Tabij)["cdik"];
    Rabij["abij"] += 0.5*       Cij["li"] * (*Tabij)["ablj"];

    // 1st term sixth line (Multiplied by *0.5 in DCD)
    Cab["da"]      = 1.0*(*Vabij)["cdkl"] * (*Tabij)["cakl"];
    Rabij["abij"] -= 1.0*       Cab["da"] * (*Tabij)["dbij"];

    // 2nd term sixth line (Multiplied by *0.5 in DCD)
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

    // calculate the amplitdues from the residuum
    doublesAmplitudesFromResiduum(Rabij);
    // and append them to the mixer
    TabijMixer->append(Rabij);
  }
}
