#include <algorithms/CcdEnergyFromCoulombFactors.hpp>
#include <math/MathFunctions.hpp>
#include <util/DryTensor.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CcdEnergyFromCoulombFactors);

CcdEnergyFromCoulombFactors::CcdEnergyFromCoulombFactors(
  std::vector<Argument> const &argumentList
): ClusterDoublesAlgorithm(argumentList) {
}

CcdEnergyFromCoulombFactors::~CcdEnergyFromCoulombFactors() {
}

//////////////////////////////////////////////////////////////////////
// Hiarata iteration routine for the CCD amplitudes Tabij (Table. 1)
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
//////////////////////////////////////////////////////////////////////
void CcdEnergyFromCoulombFactors::iterate(int i) {
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
      // For first iteration compute only the MP2 amplitudes 
      // Since Tabij = 0, Vabij is the only non-zero term
      Rabij["abij"] += (*Vabij)["abij"];
    } 
    else {
      // For the rest iterations compute the CCD amplitudes

      // Read the Coulomb Integrals Vaibj Vijkl
      Tensor<> *Vaibj(getTensorArgument("PHPHCoulombIntegrals"));
      Tensor<> *Vijkl(getTensorArgument("HHHHCoulombIntegrals"));

      // Compute the No,Nv
      int No(Vabij->lens[2]);
      int Nv(Vabij->lens[0]);

      // Define symmetries used by intermediates
      int syms[] = { NS, NS, NS, NS };

      {
	// Define intermediates
	int voov[] = { Nv, No, No, Nv };
	int vv[] = { Nv, Nv };
	int oo[] = { No, No };

	Tensor<> Kac(2, vv, syms, *Vabij->wrld, "Kac");
	Tensor<> Kki(2, oo, syms, *Vabij->wrld, "Kki");

	Tensor<> Xklij(false, *Vijkl);
	Xklij.set_name("Xklij");
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
	Rabij["abij"]  = ( 1.0) * Kac["ac"] * (*Tabij)["cbij"];

	// Contract Kki with T2 Amplitudes
	Rabij["abij"] += (-1.0) * Kki["ki"] * (*Tabij)["abkj"];

	// Build Xakic
	Xakic["akic"]  = ( 1.0) * (*Vabij)["acik"];
	Xakic["akic"] += (-0.5) * (*Vabij)["dclk"] * (*Tabij)["dail"];
	Xakic["akic"] += ( 1.0) * (*Vabij)["dclk"] * (*Tabij)["adil"];
	Xakic["akic"] += (-0.5) * (*Vabij)["cdlk"] * (*Tabij)["adil"];

	// Build Xakci
	Xakci["akci"]  = ( 1.0) * (*Vaibj)["akci"];
	Xakci["akci"] += (-0.5) * (*Vabij)["cdlk"] * (*Tabij)["dail"];

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
	Xklij["klij"]  = (*Vijkl)["klij"];
	Xklij["klij"] += (*Vabij)["cdkl"] * (*Tabij)["cdij"];

	// Contract Xklij with T2 Amplitudes
	Rabij["abij"] += Xklij["klij"] * (*Tabij)["abkl"];
      }
      
      {
	// Contract Vabcd with T2 Amplitudes via Coulomb factors

	// Read the Coulomb Factors PiqR and LambdaGR
	Tensor<complex> *PiqR(getTensorArgument<complex>("FactorOrbitals"));
	PiqR->set_name("PiqR");
	Tensor<complex> *LambdaGR(getTensorArgument<complex>("CoulombFactors"));
	LambdaGR->set_name("LambdaGR");

	int Np=No+Nv;
	int NR(PiqR->lens[1]);
	int Rvoo[] = { NR, Nv, No, No };
	int RRoo[] = { NR, NR, No, No };
	int RR[] = { NR, NR };

	Tensor<complex> VRS(2, RR, syms, *Vabij->wrld, "VRS");

	Tensor<> realXRaij(4, Rvoo, syms, *Vabij->wrld, "RealXRaij");
	Tensor<> imagXRaij(4, Rvoo, syms, *Vabij->wrld, "ImagXRaij");

	// Allocate and compute PiaR
	int aRStart[] = {No , 0};
	int aREnd[]   = {Np ,NR};
	Tensor<complex> PiaR(PiqR->slice(aRStart,aREnd));
	PiaR.set_name("PiaR");

	// Split PiaR into real and imaginary parts
	Tensor<> realPiaR(2, PiaR.lens, PiaR.sym, *PiaR.wrld, "RealPiaR");
	Tensor<> imagPiaR(2, PiaR.lens, PiaR.sym, *PiaR.wrld, "ImagPiaR");
	fromComplexTensor(PiaR, realPiaR, imagPiaR);

	// FIXME: Currently assuming GammaGqr = PiqR*PirR*LambdaGR
	//        First Pi not conjugated.
	realXRaij["Rdij"] = +1.0 * (*Tabij)["cdij"] * realPiaR["cR"];
	imagXRaij["Rdij"] = -1.0 * (*Tabij)["cdij"] * imagPiaR["cR"];
	Tensor<complex> XRaij(4, Rvoo, syms, *Vabij->wrld, "XRaij");
	toComplexTensor(realXRaij, imagXRaij, XRaij);

	Tensor<complex> XRSij(4, RRoo, syms, *Vabij->wrld, "XRSij");
	XRSij["RSij"] = XRaij["Rdij"] * PiaR["dS"];

	Univar_Function<complex> fConj(&cc4s::conj<complex>);
	Tensor<complex> conjLambdaGR(false, *LambdaGR);
	// conjLambdaGR["GR"] = conj(LambdaGR["GR"])
	conjLambdaGR.set_name("ConjLambdaGR");
	conjLambdaGR.sum(1.0, *LambdaGR,"GR", 0.0,"GR", fConj);
	VRS["RS"] = conjLambdaGR["GR"] * (*LambdaGR)["GS"];

	XRSij["RSij"] = XRSij["RSij"] * VRS["RS"];
	XRaij["Rbij"] = XRSij["RSij"]  * PiaR["bS"];

	fromComplexTensor(XRaij, realXRaij, imagXRaij);
	Rabij["abij"] += realXRaij["Rbij"]  * realPiaR["aR"];
	Rabij["abij"] += imagXRaij["Rbij"]  * imagPiaR["aR"];
      }

    }

    // Calculate the amplitdues from the residuum
    doublesAmplitudesFromResiduum(Rabij);
    // And append them to the mixer
    TabijMixer->append(Rabij);
  }
}

void CcdEnergyFromCoulombFactors::dryIterate() {
  {
    // TODO: the Mixer should provide a DryTensor in the future
    // Read the CCD amplitudes Tabij
    // DryTensor<> *Tabij(
    getTensorArgument<double, DryTensor<double>>("CcdDoublesAmplitudes");
    // );

    // Read the Coulomb Integrals Vabij Vaibj Vijkl
    DryTensor<> *Vabij(getTensorArgument<double, DryTensor<double>>("PPHHCoulombIntegrals"));
    DryTensor<> *Vaibj(getTensorArgument<double, DryTensor<double>>("PHPHCoulombIntegrals"));
    DryTensor<> *Vijkl(getTensorArgument<double, DryTensor<double>>("HHHHCoulombIntegrals"));

    // Read the Coulomb Factors PiqR and LambdaGR
    DryTensor<complex> *PiqR(getTensorArgument<complex, 
			     DryTensor<complex>>("FactorOrbitals"));
    DryTensor<complex> *LambdaGR(getTensorArgument<complex,
				 DryTensor<complex>>("CoulombFactors"));

    // Read the Particle/Hole Eigenenergies epsi epsa
    DryTensor<> *epsi(getTensorArgument<double, DryTensor<double>>("HoleEigenEnergies"));
    DryTensor<> *epsa(getTensorArgument<double, DryTensor<double>>("ParticleEigenEnergies"));
  
    // Compute the No,Nv,NR
    int No(epsi->lens[0]);
    int Nv(epsa->lens[0]);
    int NR(PiqR->lens[1]);

    int syms[] = { NS, NS, NS, NS };
    int voov[] = { Nv, No, No, Nv };
    int Rvoo[] = { NR, Nv, No, No };
    int RRoo[] = { NR, NR, No, No };
    int RR[] = { NR, NR };
    int vv[] = { Nv, Nv };
    int vR[] = { Nv, NR };
    int oo[] = { No, No };

    // Allocate Tensors for T2 amplitudes
    DryTensor<> Rabij(*Vabij);

    // Define intermediates
    DryTensor<> Kac(2, vv, syms);
    DryTensor<> Kki(2, oo, syms);

    DryTensor<> Xklij(*Vijkl);
    DryTensor<> Xakci(*Vaibj);
    DryTensor<> Xakic(4, voov, syms);

    // Contract Vabcd with T2 Amplitudes
    {
      DryTensor<complex> VRS(2, RR, syms);

      DryTensor<> realXRaij(4, Rvoo, syms);
      DryTensor<> imagXRaij(4, Rvoo, syms);

      // Allocate PiaR
      DryTensor<complex> PiaR(2, vR, syms);

      // Split PiaR into real and imaginary parts
      DryTensor<> realPiaR(2, vR, syms);
      DryTensor<> imagPiaR(2, vR, syms);

      DryTensor<complex> XRaij(4, Rvoo, syms);

      DryTensor<complex> XRSij(4, RRoo, syms);

      DryTensor<complex> conjLambdaGR(*LambdaGR);
    }

    dryDoublesAmplitudesFromResiduum(Rabij);
  }
}

