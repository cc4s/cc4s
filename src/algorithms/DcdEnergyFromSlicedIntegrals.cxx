#include <algorithms/DcdEnergyFromSlicedIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <util/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DcdEnergyFromSlicedIntegrals);

DcdEnergyFromSlicedIntegrals::DcdEnergyFromSlicedIntegrals(
  std::vector<Argument> const &argumentList
): ClusterDoublesAlgorithm(argumentList) {
}

DcdEnergyFromSlicedIntegrals::~DcdEnergyFromSlicedIntegrals() {
}

//////////////////////////////////////////////////////////////////////
// Hiarata iteration routine for the DCD amplitudes Tabij (Table. 1)
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
// Modified according to D. Kats, J. Chem. Phys. 139, 021102 (2013)
//////////////////////////////////////////////////////////////////////
void DcdEnergyFromSlicedIntegrals::iterate(int i) {
  {
    // Read the DCD amplitudes Tabij
    Tensor<> *Tabij(&TabijMixer->getNext());
    Tabij->set_name("Tabij");

    // Read the Coulomb Integrals Vabij
    Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));

    // Allocate Tensor for T2 amplitudes
    Tensor<> Rabij(false, *Vabij);
    Rabij.set_name("Rabij");

    std::string abbreviation(getAbbreviation());
    std::transform(abbreviation.begin(), abbreviation.end(), 
		   abbreviation.begin(), ::toupper);

    LOG(1, abbreviation) << "Solving T2 Amplitude Equations" << std::endl;

    if (i == 0) {
      // For first iteration compute only the MP2 amplitudes 
      // Since Tabij = 0, Vabij is the only non-zero term
      Rabij["abij"] += (*Vabij)["abij"];
    } 
    else {
      // For the rest iterations compute the DCD amplitudes

      // Read the Coulomb Integrals Vabcd Vaibj Vijkl
      Tensor<> *Vaibj(getTensorArgument("PHPHCoulombIntegrals"));
      Tensor<> *Vijkl(getTensorArgument("HHHHCoulombIntegrals"));

      // Compute the No,Nv
      int No(Vabij->lens[2]);
      int Nv(Vabij->lens[0]);

      {
	// Define intermediates
	int syms[] = { NS, NS, NS, NS };
	int voov[] = { Nv, No, No, Nv };
	int vv[] = { Nv, Nv };
	int oo[] = { No, No };

	Tensor<> Kac(2, vv, syms, *Vabij->wrld, "Kac");
	Tensor<> Kki(2, oo, syms, *Vabij->wrld, "Kki");

	//	Tensor<> Xklij(false, *Vijkl); // Removed in DCD
	//	Xklij.set_name("Xklij");       // Removed in DCD
	Tensor<> Xakci(false, *Vaibj);
	Xakci.set_name("Xakci");
	Tensor<> Xakic(4, voov, syms, *Vabij->wrld, "Xakic");

	// Build Kac
	Kac["ac"]  = (-2.0) * (*Vabij)["cdkl"] * (*Tabij)["adkl"];
	Kac["ac"] += ( 1.0) * (*Vabij)["dckl"] * (*Tabij)["adkl"];

	// Build Kki
	Kki["ki"]  = ( 2.0) * (*Vabij)["cdkl"] * (*Tabij)["cdil"];
	Kki["ki"] += (-1.0) * (*Vabij)["dckl"] * (*Tabij)["cdil"];
    
	// Contract Kac with T2 Amplitudes
	Rabij["abij"]  = ( 0.5) * Kac["ac"] * (*Tabij)["cbij"]; // Multiplied by 0.5 in DCD

	// Contract Kki with T2 Amplitudes
	Rabij["abij"] += (-0.5) * Kki["ki"] * (*Tabij)["abkj"]; // Multiplied by 0.5 in DCD

	// Build Xakic
	Xakic["akic"]  = ( 1.0) * (*Vabij)["acik"];
	Xakic["akic"] += (-0.5) * (*Vabij)["dclk"] * (*Tabij)["dail"];
	Xakic["akic"] += ( 1.0) * (*Vabij)["dclk"] * (*Tabij)["adil"];
	//	Xakic["akic"] += (-0.5) * (*Vabij)["cdlk"] * (*Tabij)["adil"]; // Removed in DCD

	// Build Xakci
	Xakci["akci"]  = ( 1.0) * (*Vaibj)["akci"];
	//	Xakci["akci"] += (-0.5) * (*Vabij)["cdlk"] * (*Tabij)["dail"]; // Removed in DCD

	// Contract Xakic and Xakci intermediates with T2 amplitudes Tabij
	Rabij["abij"] += ( 2.0) * Xakic["akic"] * (*Tabij)["cbkj"];
	Rabij["abij"] += (-1.0) * Xakic["akic"] * (*Tabij)["bckj"];

	Rabij["abij"] += (-1.0) * Xakci["akci"] * (*Tabij)["cbkj"];
	Rabij["abij"] += (-1.0) * Xakci["bkci"] * (*Tabij)["ackj"];

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
	//	Xklij["klij"]  = (*Vijkl)["klij"]; // Removed in DCD
	//	Xklij["klij"] += (*Vabij)["cdkl"] * (*Tabij)["cdij"]; // Removed in DCD

	// Contract Xklij with T2 Amplitudes
	Rabij["abij"] += (*Vijkl)["klij"] * (*Tabij)["abkl"]; // Xklij replaced by Vklij in DCD
      }
      
      // Contract Vabcd with T2 Amplitudes
      {
	// Slice Vabcd
	// Read the sliceRank. If not provided use No
	int sliceRank(getIntegerArgument
		      ("sliceRank",No));

	// Slice loop starts here
	for (int b(0); b < Nv; b += sliceRank) {
	  for (int a(b); a < Nv; a += sliceRank) {
	    LOG(1, abbreviation) << "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
	    Tensor<> *Vxycd(sliceCoulombIntegrals(a, b, sliceRank));
	    Vxycd->set_name("Vxycd");
	    int lens[] = { Vxycd->lens[0], Vxycd->lens[1], No, No };
	    int syms[] = {NS, NS, NS, NS};
	    Tensor<> Rxyij(4, lens, syms, *Vxycd->wrld, "Rxyij");
	    Rxyij["xyij"] = (*Vxycd)["xycd"] * (*Tabij)["cdij"];
	    sliceIntoResiduum(Rxyij, a, b, Rabij);
	    // The integrals of this slice are not needed anymore
	    delete Vxycd;
	  }
	}

      }

    }

    // Calculate the amplitdues from the residuum
    doublesAmplitudesFromResiduum(Rabij);
    // And append them to the mixer
    TabijMixer->append(Rabij);
  }
}

void DcdEnergyFromSlicedIntegrals::dryIterate() {
  {
    // TODO: the Mixer should provide a DryTensor in the future
    // Read the DCD amplitudes Tabij
    getTensorArgument<double, DryTensor<double>>("DcdDoublesAmplitudes");

    // Read the Coulomb Integrals Vabcd Vabij Vaibj Vijkl
    DryTensor<> *Vabij(getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals"));
    DryTensor<> *Vaibj(getTensorArgument<double, DryTensor<double>>("PHPHCoulombIntegrals"));
    DryTensor<> *Vijkl(getTensorArgument<double, DryTensor<double>>("HHHHCoulombIntegrals"));
  
    // Compute the No,Nv,Np
    int No(Vabij->lens[2]);
    int Nv(Vabij->lens[0]);

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

    {
      // Read the sliceRank. If not provided use No
      int sliceRank(getIntegerArgument
		    ("sliceRank",No));

      DryTensor<> *Vxycd(drySliceCoulombIntegrals(sliceRank));
      int lens[] = { Vxycd->lens[0], Vxycd->lens[1], No, No };
      int syms[] = {NS, NS, NS, NS};
      DryTensor<> Rxyij(4, lens, syms);
    }

    dryDoublesAmplitudesFromResiduum(Rabij);
  }
}

//////////////////////////////////////////////////////////////////////
// Bartlett iteration routine for the DCD amplitudes Tabij 
// Rev. Mod. Phys. 79, 291  Page 305, Figure 8. -> CCD
// J. Chem. Phys. 139, 021102 (2013) -> DCD
//////////////////////////////////////////////////////////////////////
void DcdEnergyFromSlicedIntegrals::iterateBartlett(int i) {
  {
    // Read the CCD amplitudes Tabij
    Tensor<> *Tabij(&TabijMixer->getNext());
    Tabij->set_name("Tabij");

    // Read the Coulomb Integrals Vabij
    Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));

    // Allocate Tensor for T2 amplitudes
    Tensor<> Rabij(false, *Vabij);
    Rabij.set_name("Rabij");

    std::string abbreviation(getAbbreviation());
    std::transform(abbreviation.begin(), abbreviation.end(), 
		   abbreviation.begin(), ::toupper);
  
    LOG(1, abbreviation) << "Solving T2 Amplitude Equations" << std::endl;

    if (i == 0) {
      // For first iteration compute only the MP2 amplitudes since Tabij = 0

      // Rabij contracted with Vabij is the only non-zero term
      Rabij["abij"] += (*Vabij)["abij"];
      
    } 
    else {
      // For the rest iterations compute the CCD amplitudes

      // Read the Coulomb Integrals Vabcd Vabij Vaibj Vijkl
      Tensor<> *Vaibj(getTensorArgument("PHPHCoulombIntegrals"));
      Tensor<> *Vijkl(getTensorArgument("HHHHCoulombIntegrals"));

      // Compute the No,Nv
      int No(Vabij->lens[2]);
      int Nv(Vabij->lens[0]);

      {
	//////////////////////////////////////////////////////////////////////
	// Create linear terms with T2 Amplitudes that need permutation
	//////////////////////////////////////////////////////////////////////

	// Contract Vabcd with T2 Amplitudes (4th term first line)
	{
	  // Slice Vabcd
	  // Read the sliceRank. If not provided use No
	  int64_t sliceRank(getIntegerArgument
			    ("sliceRank",No));

	  // Slice loop starts here
	  for (int b(0); b < Nv; b += sliceRank) {
	    for (int a(b); a < Nv; a += sliceRank) {
	      LOG(1, abbreviation) << "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
	      Tensor<> *Vxyef(sliceCoulombIntegrals(a, b, sliceRank));
	      int lens[] = { Vxyef->lens[0], Vxyef->lens[1], No, No };
	      int syms[] = {NS, NS, NS, NS};
	      Tensor<> Rxyij(4, lens, syms, *Vxyef->wrld, "Rxyij");
	      Rxyij["xyij"] = ( 0.5) * (*Vxyef)["xyef"] * (*Tabij)["efij"];
	      sliceIntoResiduum(Rxyij, a, b, Rabij);
	      // The integrals of this slice are not needed anymore
	      delete Vxyef;
	    }
	  }

	}

	// Contract Vijkl with T2 Amplitudes (5th term first line)
	Rabij["abij"] += ( 0.5) * (*Vijkl)["mnij"] * (*Tabij)["abmn"];

	// Contract Vabij with T2 Amplitudes (1st term second line)
	Rabij["abij"] += ( 2.0) * (*Vabij)["ebmj"] * (*Tabij)["aeim"];

	// Contract Vabij with T2 Amplitudes (2nd term second line)
	Rabij["abij"] += (-1.0) * (*Vabij)["ebmj"] * (*Tabij)["eaim"];

	// Contract Vaibj with T2 Amplitudes (3rd term second line)
	Rabij["abij"] += (-1.0) * (*Vaibj)["eibm"] * (*Tabij)["aemj"];

	// Contract Vaibj with T2 Amplitudes (4th term second line)
	Rabij["abij"] += (-1.0) * (*Vaibj)["ejbm"] * (*Tabij)["aeim"];

	//////////////////////////////////////////////////////////////////////
	// Create quadratic terms with T2 Amplitudes that need permutation
	//////////////////////////////////////////////////////////////////////

	// 1st term third line
	Rabij["abij"] += ( 2.0) * (*Tabij)["aeim"] * (*Vabij)["efmn"] * (*Tabij)["fbnj"];

	// 2nd term third line
	Rabij["abij"] += (-2.0) * (*Tabij)["aeim"] * (*Vabij)["efmn"] * (*Tabij)["fbjn"];

	// 3rd term third line
	Rabij["abij"] += ( 0.5) * (*Tabij)["eaim"] * (*Vabij)["efmn"] * (*Tabij)["fbjn"];

	// 1st term fourth line
	//	Rabij["abij"] += (-1.0) * (*Tabij)["aeim"] * (*Vabij)["femn"] * (*Tabij)["fbnj"]; // Removed in DCD

	// 2nd term fourth line
	//	Rabij["abij"] += ( 1.0) * (*Tabij)["aemi"] * (*Vabij)["femn"] * (*Tabij)["fbnj"]; // Removed in DCD

	// 3rd term fourth line
	//	Rabij["abij"] += ( 0.5) * (*Tabij)["aemj"] * (*Vabij)["femn"] * (*Tabij)["fbin"]; // Removed in DCD

	// 1st term fifth line
	//	Rabij["abij"] += ( 0.5) * (*Tabij)["abmn"] * (*Vabij)["efmn"] * (*Tabij)["efij"]; // Removed in DCD

	// 2nd term fifth line (Multiplied by *0.5 in DCD)
	Rabij["abij"] += (-1.0) * (*Tabij)["abnj"] * (*Vabij)["efmn"] * (*Tabij)["efmi"];

	// 3rd term fifth line (Multiplied by *0.5 in DCD)
	Rabij["abij"] += ( 0.5) * (*Tabij)["abnj"] * (*Vabij)["efmn"] * (*Tabij)["efim"];

	// 1st term sixth line (Multiplied by *0.5 in DCD)
	Rabij["abij"] += (-1.0) * (*Tabij)["fbij"] * (*Vabij)["efmn"] * (*Tabij)["eamn"];

	// 2nd term sixth line (Multiplied by *0.5 in DCD)
	Rabij["abij"] += ( 0.5) * (*Tabij)["fbij"] * (*Vabij)["efmn"] * (*Tabij)["aemn"];
      }

      //////////////////////////////////////////////////////////////////////
      // Symmetrize Rabij by applying permutation operator
      // To save memory we use Caibj as intermediate for the permutation operator 
      //////////////////////////////////////////////////////////////////////
    
      {
	// Tensor used for permutation operation
	Tensor<> Caibj(false, *Vaibj);

	Caibj["aibj"]  = Rabij["abij"];
	Rabij["abij"] += Caibj["bjai"];
      }

      //////////////////////////////////////////////////////////////////////
      // Add the Vabij (the only term that does not need permutation)
      //////////////////////////////////////////////////////////////////////
      Rabij["abij"] += (*Vabij)["abij"];
    }

    // calculate the amplitdues from the residuum
    doublesAmplitudesFromResiduum(Rabij);
    // and append them to the mixer
    TabijMixer->append(Rabij);
  }
}