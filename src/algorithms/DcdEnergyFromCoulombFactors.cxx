#include <algorithms/DcdEnergyFromCoulombFactors.hpp>
#include <math/MathFunctions.hpp>
#include <math/Complex.hpp>
#include <math/ComplexTensor.hpp>
#include <mixers/Mixer.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <Cc4s.hpp>
#include <ctf.hpp>
#include <Options.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DcdEnergyFromCoulombFactors);

DcdEnergyFromCoulombFactors::DcdEnergyFromCoulombFactors(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

DcdEnergyFromCoulombFactors::~DcdEnergyFromCoulombFactors() {
  // delete the mixer to free its resources
  if (TabijMixer) delete TabijMixer;
}

/**
 * \brief Calculates DCD energy from Coulomb integrals Vabcd Vabij Vaibj Vijkl
 */

void DcdEnergyFromCoulombFactors::run() {
  // Read the Coulomb Integrals Vabij required for the DCD energy
  Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));

  // Read the Particle/Hole Eigenenergies epsi epsa required for the DCD energy
  Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
  Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  
  // Compute the no,nv,np
  int no(epsi->lens[0]);
  int nv(epsa->lens[0]);

  // instantiate mixer for the amplitudes, by default use the trivial mixer
  std::string mixerName(getTextArgument("mixer", "TrivialMixer"));
  TabijMixer = MixerFactory<double>::create(mixerName);
  if (!TabijMixer) {
    std::stringstream stringStream;
    stringStream << "Mixer not implemented: " << mixerName;
    throw new Exception(stringStream.str());
  }
  {
    // Allocate the DCD amplitudes Tabij and append it to the mixer
    int syms[] = { NS, NS, NS, NS };
    int vvoo[] = { nv, nv, no, no };
    Tensor<> Tabij(4, vvoo, syms, *Cc4s::world, "Tabij");
    TabijMixer->append(Tabij);
    // the amplitudes will from now on be managed by the mixer
  }

  // Allocate the DCD energy e
  Scalar<> energy(*Cc4s::world);
  double e(0), dire, exce;

  LOG(0, "DCD") <<
    "Solving Distinguishable Cluster Doubles Amplitude Equations:" <<
    std::endl;

  // Iteration for determining the DCD amplitudes Tabij
  // and the Dcd energy e
  int64_t maxIterationsCount(
    getIntegerArgument("maxIterations", DEFAULT_MAX_ITERATIONS)
  );
  for (int i(0); i < maxIterationsCount; ++i) {
    LOG(0, "DCD") << "iteration: " << i+1 << std::endl;
    iterate(i);
    // get the current amplitdues from the mixer
    Tensor<> *Tabij(&TabijMixer->getNext());
    energy[""] = 2.0 * (*Tabij)["abij"] * (*Vabij)["abij"];
    dire = energy.get_val();
    energy[""] = (*Tabij)["abji"] * (*Vabij)["abij"];
    exce = -1.0 * energy.get_val();
    e = dire + exce;
    LOG(0, "DCD") << "e=" << e << std::endl;
    LOG(1, "DCD") << "DCDdir=" << dire << std::endl;
    LOG(1, "DCD") << "DCDexc=" << exce << std::endl;
  }

  LOG(1, "DCD") << "DCD correlation energy = " << e << std::endl;

  // return a copy of the amplitudes contained in the mixer
  // the amplitdues contained in the mixer will be deleted upon
  // destruction of this algorithm
  allocatedTensorArgument(
    "DcdDoublesAmplitudes", new Tensor<>(TabijMixer->getNext())
  );
  setRealArgument("DcdEnergy", e);
}

//////////////////////////////////////////////////////////////////////
// Hiarata iteration routine for the DCD amplitudes Tabij (Table. 1)
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
// Modified according to D. Kats, J. Chem. Phys. 139, 021102 (2013)
//////////////////////////////////////////////////////////////////////
void DcdEnergyFromCoulombFactors::iterate(int i) {
  {
    // get the current amplitdues for this iteration from the mixer
    Tensor<> *Tabij(&TabijMixer->getNext());

    // Read the Coulomb Integrals Vabcd Vabij Vaibj Vijkl
    Tensor<> *Vabij(getTensorArgument("PPHHCoulombIntegrals"));
    Tensor<> *Vaibj(getTensorArgument("PHPHCoulombIntegrals"));
    Tensor<> *Vijkl(getTensorArgument("HHHHCoulombIntegrals"));

    // Read the Coulomb Factors PiqR and LambdaGR
    Tensor<complex> *PiqR(getTensorArgument<complex>("FactorOrbitals"));
    Tensor<complex> *LambdaGR(getTensorArgument<complex>("CoulombFactors"));

    // Read the Particle/Hole Eigenenergies epsi epsa
    Tensor<> *epsi(getTensorArgument<>("HoleEigenEnergies"));
    Tensor<> *epsa(getTensorArgument<>("ParticleEigenEnergies"));
  
    // Compute the no,nv,np
    int no(epsi->lens[0]);
    int nv(epsa->lens[0]);
    int np=no+nv;
    int nR(PiqR->lens[1]);

    int syms[] = { NS, NS, NS, NS };
    int voov[] = { nv, no, no, nv };
    int Rvoo[] = { nR, nv, no, no };
    int RRoo[] = { nR, nR, no, no };
    int RR[] = { nR, nR };
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
    // Rabij["abij"] += (*Vabcd)["abcd"] * (*Tabij)["cdij"];
    {
      Tensor<complex> VRS(2, RR, syms, *Cc4s::world, "VRS");

      Tensor<> realXRaij(4, Rvoo, syms, *Cc4s::world, "RealXRaij");
      Tensor<> imagXRaij(4, Rvoo, syms, *Cc4s::world, "ImagXRaij");

      // Allocate and compute PiaR
      int aRStart[] = {no , 0};
      int aREnd[]   = {np ,nR};
      Tensor<complex> PiaR(PiqR->slice(aRStart,aREnd));

      // Split PiaR into real and imaginary parts
      Tensor<> realPiaR(
        2, PiaR.lens, PiaR.sym, *PiaR.wrld, "RealPiaR"
      );
      Tensor<> imagPiaR(
        2, PiaR.lens, PiaR.sym, *PiaR.wrld, "ImagPiaR"
      );
      fromComplexTensor(PiaR, realPiaR, imagPiaR);

      // FIXME: Currently assuming GammaGqr = PiqR*PirR*LambdaGR
      //        First Pi not conjugated.
      realXRaij["Rdij"] = +1.0 * (*Tabij)["cdij"] * realPiaR["cR"];
      imagXRaij["Rdij"] = -1.0 * (*Tabij)["cdij"] * imagPiaR["cR"];
      Tensor<complex> XRaij(4, Rvoo, syms, *Cc4s::world, "XRaij");
      toComplexTensor(realXRaij, imagXRaij, XRaij);

      Tensor<complex> XRSij(4, RRoo, syms, *Cc4s::world, "XRSij");
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

    // Build Dabij
    Dabij["abij"]  = (*epsi)["i"];
    Dabij["abij"] += (*epsi)["j"];
    Dabij["abij"] -= (*epsa)["a"];
    Dabij["abij"] -= (*epsa)["b"];
    // NOTE: ctf double counts if lhs tensor is SH,SH
    Dabij["abij"] = Dabij["abij"];

    // Divide Rabij/Dabij to get Tabij
    Bivar_Function<> fDivide(&divide<double>);
    // dont't use Tabij, since it actually is part of the mixer
    // and must be modified. Rabij is, however, available
    Rabij.contract(1.0, Rabij,"abij", Dabij,"abij", 0.0,"abij", fDivide);

    // pass the new amplitudes to the mixer
    TabijMixer->append(Rabij);
  }
}

