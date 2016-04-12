#include <algorithms/DcdEnergyFromCoulombFactors.hpp>
#include <math/MathFunctions.hpp>
#include <util/DryTensor.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DcdEnergyFromCoulombFactors);

DcdEnergyFromCoulombFactors::DcdEnergyFromCoulombFactors(
  std::vector<Argument> const &argumentList
): ClusterDoublesAlgorithm(argumentList) {
}

DcdEnergyFromCoulombFactors::~DcdEnergyFromCoulombFactors() {
}

//////////////////////////////////////////////////////////////////////
// Hiarata iteration routine for the DCD amplitudes Tabij (Table. 1)
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
// Modified according to D. Kats, J. Chem. Phys. 139, 021102 (2013)
//////////////////////////////////////////////////////////////////////
void DcdEnergyFromCoulombFactors::iterate(int i) {
  {
    // Read the DCD amplitudes Tabij
    Tensor<> *Tabij(&TabijMixer->getNext());

    // Read the Coulomb Integrals Vabij Vaibj Vijkl
    Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
    Tensor<> *Vaibj(getTensorArgument("PHPHCoulombIntegrals"));
    Tensor<> *Vijkl(getTensorArgument("HHHHCoulombIntegrals"));

    // Read the Coulomb Factors PiqR and LambdaGR
    Tensor<complex> *PiqR(getTensorArgument<complex>("FactorOrbitals"));
    Tensor<complex> *LambdaGR(getTensorArgument<complex>("CoulombFactors"));

    // Read the Particle/Hole Eigenenergies epsi epsa
    Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
    Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  
    // Compute the No,Nv,Np,NR
    int No(epsi->lens[0]);
    int Nv(epsa->lens[0]);
    int Np=No+Nv;
    int NR(PiqR->lens[1]);

    int syms[] = { NS, NS, NS, NS };
    int voov[] = { Nv, No, No, Nv };
    int Rvoo[] = { NR, Nv, No, No };
    int RRoo[] = { NR, NR, No, No };
    int RR[] = { NR, NR };
    int vv[] = { Nv, Nv };
    int oo[] = { No, No };

    // Allocate Tensors for T2 amplitudes
    Tensor<> Rabij(false, *Vabij);

    // Define intermediates
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
    // Xklij["klij"] += (*Vabij)["cdkl"] * (*Tabij)["cdij"]; //Removed in Dcd

    // Contract Xklij with T2 Amplitudes
    Rabij["abij"] += Xklij["klij"] * (*Tabij)["abkl"];

    // Contract Vabcd with T2 Amplitudes
    {
      Tensor<complex> VRS(2, RR, syms, *epsi->wrld, "VRS");

      Tensor<> realXRaij(4, Rvoo, syms, *epsi->wrld, "RealXRaij");
      Tensor<> imagXRaij(4, Rvoo, syms, *epsi->wrld, "ImagXRaij");

      // Allocate and compute PiaR
      int aRStart[] = {No , 0};
      int aREnd[]   = {Np ,NR};
      Tensor<complex> PiaR(PiqR->slice(aRStart,aREnd));

      // Split PiaR into real and imaginary parts
      Tensor<> realPiaR(2, PiaR.lens, PiaR.sym, *PiaR.wrld, "RealPiaR");
      Tensor<> imagPiaR(2, PiaR.lens, PiaR.sym, *PiaR.wrld, "ImagPiaR");
      fromComplexTensor(PiaR, realPiaR, imagPiaR);

      // FIXME: Currently assuming GammaGqr = PiqR*PirR*LambdaGR
      //        First Pi not conjugated.
      realXRaij["Rdij"] = +1.0 * (*Tabij)["cdij"] * realPiaR["cR"];
      imagXRaij["Rdij"] = -1.0 * (*Tabij)["cdij"] * imagPiaR["cR"];
      Tensor<complex> XRaij(4, Rvoo, syms, *epsi->wrld, "XRaij");
      toComplexTensor(realXRaij, imagXRaij, XRaij);

      Tensor<complex> XRSij(4, RRoo, syms, *epsi->wrld, "XRSij");
      XRSij["RSij"] = XRaij["Rdij"] * PiaR["dS"];

      Univar_Function<complex> fConj(&cc4s::conj<complex>);
      Tensor<complex> conjLambdaGR(false, *LambdaGR);
      // conjLambdaGR["GR"] = conj(LambdaGR["GR"])
      conjLambdaGR.sum(1.0, *LambdaGR,"GR", 0.0,"GR", fConj);
      VRS["RS"] = conjLambdaGR["GR"] * (*LambdaGR)["GS"];

      XRSij["RSij"] = XRSij["RSij"] * VRS["RS"];
      XRaij["Rbij"] = XRSij["RSij"]  * PiaR["bS"];

      fromComplexTensor(XRaij, realXRaij, imagXRaij);
      Rabij["abij"] += realXRaij["Rbij"]  * realPiaR["aR"];
      Rabij["abij"] += imagXRaij["Rbij"]  * imagPiaR["aR"];
    }

    // calculate the amplitdues from the residuum
    doublesAmplitudesFromResiduum(Rabij);
    // and append them to the mixer
    TabijMixer->append(Rabij);
  }
}

void DcdEnergyFromCoulombFactors::dryIterate() {
  {
    // TODO: the Mixer should provide a DryTensor in the future
    // Read the DCD amplitudes Tabij
    // DryTensor<> *Tabij(
    getTensorArgument<double, DryTensor<double>>("DcdDoublesAmplitudes");
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
  
    // Compute the No,Nv,Np,NR
    int No(epsi->lens[0]);
    int Nv(epsa->lens[0]);
    int Np=No+Nv;
    int NR(PiqR->lens[1]);

    int syms[] = { NS, NS, NS, NS };
    int voov[] = { Nv, No, No, Nv };
    int Rvoo[] = { NR, Nv, No, No };
    int RRoo[] = { NR, NR, No, No };
    int RR[] = { NR, NR };
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

    // TODO: implment dryDoublesAmplitudesFromResiduum
    // at the moment, assume usage of Dabij
    DryTensor<> Dabij(*Vabij);
  }
}

