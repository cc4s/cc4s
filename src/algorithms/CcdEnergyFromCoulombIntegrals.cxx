#include <algorithms/CcdEnergyFromCoulombIntegrals.hpp>
#include <math/MathFunctions.hpp>
#include <tcc/DryTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CcdEnergyFromCoulombIntegrals);

CcdEnergyFromCoulombIntegrals::CcdEnergyFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): ClusterDoublesAlgorithm(argumentList) {
}

CcdEnergyFromCoulombIntegrals::~CcdEnergyFromCoulombIntegrals() {
}

//////////////////////////////////////////////////////////////////////
// Hiarata iteration routine for the CCD amplitudes Tabij (Table. 1)
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
//////////////////////////////////////////////////////////////////////
void CcdEnergyFromCoulombIntegrals::iterate(int i) {
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

    if (i == 0 && !isArgumentGiven("startingDoublesAmplitudes") ) {
      // For first iteration compute only the MP2 amplitudes 
      // Since Tabij = 0, Vabij is the only non-zero term
      Rabij["abij"] += (*Vabij)["abij"];
    } 
    else {
      // For the rest iterations compute the CCD amplitudes

      // Read the Coulomb Integrals Vabcd Vaibj Vijkl
      // the PPPPCoulombIntegrals may not be given then slicing is required
      Tensor<> *Vabcd(isArgumentGiven("PPPPCoulombIntegrals") ?
                      getTensorArgument("PPPPCoulombIntegrals") : nullptr);
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
      
      // Contract Vabcd with T2 Amplitudes
      if (Vabcd) {
        Rabij["abij"] += (*Vabcd)["abcd"] * (*Tabij)["cdij"];
      } 
      else {
        if (isArgumentGiven("CoulombFactors")) {
          // Read the factorsSliceSize. If not provided use NG.
          Tensor<complex> *LambdaGR(getTensorArgument<complex>("CoulombFactors"));
          LambdaGR->set_name("LambdaGR");

          int NR(LambdaGR->lens[1]);

          int factorsSliceSize(getIntegerArgument
                           ("factorsSliceSize",NR));


          // Allocate Tensor for T2 amplitudes
          Tensor<> Sabij(false, *Vabij);
          Sabij.set_name("Sabij");
    
          // Slice loop starts here
            for (int a(0); a < NR; a += factorsSliceSize) {
              LOG(1, abbreviation) << "Evaluting residuum from coulomb factors at R=" 
                                   << a << std::endl;
              Tensor<> *Fabij(sliceAmplitudesFromCoulombFactorsTcc(a, factorsSliceSize));
              Fabij->set_name("Fabij");
              Sabij["abij"] += (*Fabij)["abij"];
              delete Fabij;
            }
          Rabij["abij"] += 0.5 * Sabij["abij"];
          Rabij["abij"] += 0.5 * Sabij["baji"];
        }

        else {
          // Read the integralsSliceSize. If not provided use No
          int integralsSliceSize(getIntegerArgument
                        ("integralsSliceSize",No));

          // Slice loop starts here
          for (int b(0); b < Nv; b += integralsSliceSize) {
            for (int a(b); a < Nv; a += integralsSliceSize) {
              LOG(1, abbreviation) << "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
              Tensor<> *Vxycd(sliceCoulombIntegrals(a, b, integralsSliceSize));
              Vxycd->set_name("Vxycd");
              int lens[] = { Vxycd->lens[0], Vxycd->lens[1], No, No };
              int syms[] = {NS, NS, NS, NS};
              Tensor<> Rxyij(4, lens, syms, *Vxycd->wrld);
              Rxyij["xyij"] = (*Vxycd)["xycd"] * (*Tabij)["cdij"];
              sliceIntoResiduum(Rxyij, a, b, Rabij);
              // The integrals of this slice are not needed anymore
              delete Vxycd;
            }
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

void CcdEnergyFromCoulombIntegrals::dryIterate() {
  {
    // TODO: the Mixer should provide a DryTensor in the future
    // Read the CCD amplitudes Tabij
    getTensorArgument<double, DryTensor<double>>("CcdDoublesAmplitudes");

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

    std::string abbreviation(getAbbreviation());
    std::transform(abbreviation.begin(), abbreviation.end(), 
                   abbreviation.begin(), ::toupper);
  
    // Compute the No,Nv,Np
    int No(Vabij->lens[2]);
    int Nv(Vabij->lens[0]);

    int syms[] = { NS, NS, NS, NS };
    int voov[] = { Nv, No, No, Nv };
    int vv[] = { Nv, Nv };
    int oo[] = { No, No };

    // Allocate Tensors for T2 amplitudes
    DryTensor<> Rabij(*Vabij, SOURCE_LOCATION);

    // Define intermediates
    DryTensor<> Kac(2, vv, syms, SOURCE_LOCATION);
    DryTensor<> Kki(2, oo, syms, SOURCE_LOCATION);

    DryTensor<> Xklij(*Vijkl, SOURCE_LOCATION);
    DryTensor<> Xakci(*Vaibj, SOURCE_LOCATION);
    DryTensor<> Xakic(4, voov, syms, SOURCE_LOCATION);

    if (!Vabcd) {
      if (isArgumentGiven("CoulombFactors")) {
        // Read the factorsSliceSize. If not provided use NG.
        DryTensor<complex> *LambdaGR(getTensorArgument<complex, 
                               DryTensor<complex>>("CoulombFactors"));

        int NR(LambdaGR->lens[1]);

        int factorsSliceSize(getIntegerArgument("factorsSliceSize",NR));

        LOG(1, abbreviation) << "Computing residuum Rabij from factors with NR=" << NR 
                             << ", using slicing size=" << factorsSliceSize << std::endl;
        DryTensor<> *Fabij(drySliceAmplitudesFromCoulombFactors(factorsSliceSize));
        delete Fabij;
        }
        else {
          // Read the integralsSliceSize. If not provided use No
          int integralsSliceSize(getIntegerArgument
                                 ("integralsSliceSize",No));

          LOG(1, abbreviation) << "Slicing Vabcd with Nv=" << Nv << ", with integals slice size=" 
                               << integralsSliceSize << std::endl;

          // Slice if Vabcd is not specified
          DryTensor<> *Vxycd(drySliceCoulombIntegrals(integralsSliceSize));
          int lens[] = { Vxycd->lens[0], Vxycd->lens[1], No, No };
          int syms[] = {NS, NS, NS, NS};
          DryTensor<> Rxyij(4, lens, syms, SOURCE_LOCATION);
        }
    }

    dryDoublesAmplitudesFromResiduum(Rabij);
  }
}

//////////////////////////////////////////////////////////////////////
// Bartlett iteration routine for the CCD amplitudes Tabij 
// Rev. Mod. Phys. 79, 291  Page 305, Figure 8. -> CCD
//////////////////////////////////////////////////////////////////////
void CcdEnergyFromCoulombIntegrals::iterateBartlett(int i) {
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
      Tensor<> *Vabcd(getTensorArgument("PPPPCoulombIntegrals"));
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
        if (Vabcd) {
          Rabij["abij"]  = ( 0.5) * (*Vabcd)["abef"] * (*Tabij)["efij"];
        } 
        else {
          // Slice if Vabcd is not specified

          // Read the integralsSliceSize. If not provided use No
          int64_t integralsSliceSize(getIntegerArgument
                            ("integralsSliceSize",No));

          // Slice loop starts here
          for (int b(0); b < Nv; b += integralsSliceSize) {
            for (int a(b); a < Nv; a += integralsSliceSize) {
              LOG(1, abbreviation) << "Evaluting Vabcd at a=" << a << ", b=" << b << std::endl;
              Tensor<> *Vxyef(sliceCoulombIntegrals(a, b, integralsSliceSize));
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
        Rabij["abij"] += (-1.0) * (*Tabij)["aeim"] * (*Vabij)["femn"] * (*Tabij)["fbnj"];

        // 2nd term fourth line
        Rabij["abij"] += ( 1.0) * (*Tabij)["aemi"] * (*Vabij)["femn"] * (*Tabij)["fbnj"];

        // 3rd term fourth line
        Rabij["abij"] += ( 0.5) * (*Tabij)["aemj"] * (*Vabij)["femn"] * (*Tabij)["fbin"];

        // 1st term fifth line
        Rabij["abij"] += ( 0.5) * (*Tabij)["abmn"] * (*Vabij)["efmn"] * (*Tabij)["efij"];

        // 2nd term fifth line
        Rabij["abij"] += (-2.0) * (*Tabij)["abnj"] * (*Vabij)["efmn"] * (*Tabij)["efmi"];

        // 3rd term fifth line
        Rabij["abij"] += ( 1.0) * (*Tabij)["abnj"] * (*Vabij)["efmn"] * (*Tabij)["efim"];

        // 1st term sixth line
        Rabij["abij"] += (-2.0) * (*Tabij)["fbij"] * (*Vabij)["efmn"] * (*Tabij)["eamn"];

        // 2nd term sixth line
        Rabij["abij"] += ( 1.0) * (*Tabij)["fbij"] * (*Vabij)["efmn"] * (*Tabij)["aemn"];
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
