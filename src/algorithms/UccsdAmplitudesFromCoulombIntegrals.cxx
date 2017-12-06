#include <algorithms/UccsdAmplitudesFromCoulombIntegrals.hpp>
#include <unistd.h>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <math/RandomTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>

using namespace cc4s;
using namespace tcc;

ALGORITHM_REGISTRAR_DEFINITION(UccsdAmplitudesFromCoulombIntegrals);

UccsdAmplitudesFromCoulombIntegrals::UccsdAmplitudesFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): ClusterSinglesDoublesAlgorithm(argumentList) {
}


UccsdAmplitudesFromCoulombIntegrals::~UccsdAmplitudesFromCoulombIntegrals() {
}

bool UccsdAmplitudesFromCoulombIntegrals::maskIsGiven() {
  return isArgumentGiven("MaskVirtualRange") &&
          isArgumentGiven("MaskParticleRange");
}

void UccsdAmplitudesFromCoulombIntegrals::run() {
  if (maskIsGiven()) {
    LOG(0, getAbbreviation()) << "Mask selected" << std::endl;
    createMask();
  }
  ClusterSinglesDoublesAlgorithm::run();
}

void UccsdAmplitudesFromCoulombIntegrals::createMask(){
  LOG(1, getAbbreviation()) << "Creating mask" << std::endl;

  auto epsi(getTensorArgument<double>("HoleEigenEnergies"));
  auto epsa(getTensorArgument<double>("ParticleEigenEnergies"));
  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int vo[] = {Nv, No}, vvoo[] = {Nv,Nv,No,No};
  int syms[] = {NS, NS};
  int syms_4[] = {NS, NS, NS, NS};

  RangeParser virtualRange(getTextArgument("MaskVirtualRange"));
  RangeParser particleRange(getTextArgument("MaskParticleRange"));
  LOG(0, getAbbreviation()) << "Nv, No: " << Nv << ", " << No << std::endl;
  LOG(0, getAbbreviation()) << "Virt. Range: " << virtualRange;
  LOG(0, getAbbreviation()) << "Part. Range: " << particleRange;

  if (
    virtualRange.get_max()  >= Nv ||
    particleRange.get_max() >= No
  ) {
    LOG(0, getAbbreviation()) << "Mask range out of bounds" << std::endl;
    throw new EXCEPTION("Mask range out of bounds");
  }

  Mai = NEW(CTF::Tensor<double>, 2, vo, syms, *Cc4s::world, "Mai");
  Mabij = NEW(CTF::Tensor<double>, 4, vvoo, syms_4, *Cc4s::world, "Mabij");
  (*Mai)["ai"] = 1.0;
  LOG(0, getAbbreviation()) << "Mai done" << std::endl;
  (*Mabij)["abij"] = 1.0;
  LOG(0, getAbbreviation()) << "Mabij done" << std::endl;

  int64_t *MaiIndex, *MabijIndex;
  double *MaiValue, *MabijValue;
  int64_t MaiCount, MabijCount;

  LOG(1, getAbbreviation()) << "Writing Mai and Mabij" << std::endl;
  for (auto a : virtualRange.getRange()) {
  for (auto i : particleRange.getRange()) {
    if (Mai->wrld->rank == 0) {
      MaiCount = 1;
      MaiIndex = (int64_t*) malloc(MaiCount);
      MaiValue = (double*) malloc(MaiCount);
      MaiIndex[0] = 0.0;
      MaiIndex[0] = a + i*Nv;
    } else {
      MaiCount = 0;
      MaiIndex = (int64_t*) malloc(MaiCount);
      MaiValue = (double*) malloc(MaiCount);
    }
    Mai->write(MaiCount, MaiIndex, MaiValue);
  for (auto b : virtualRange.getRange()) {
  for (auto j : particleRange.getRange()) {
    if (Mabij->wrld->rank == 0) {
      MabijCount = 1;
      MabijIndex = (int64_t*) malloc(MabijCount);
      MabijValue = (double*) malloc(MabijCount);
      MabijIndex[0] = 0.0;
      MabijIndex[0] = a + b*Nv + i*Nv*Nv + j*Nv*Nv*No;
    } else {
      MabijCount = 0;
      MabijIndex = (int64_t*) malloc(MabijCount);
      MabijValue = (double*) malloc(MabijCount);
    }
    Mabij->write(MabijCount, MabijIndex, MabijValue);
  }
  }
  }
  }

  //Mai->print();
  //Mabij->print();
}

PTR(FockVector<complex>) UccsdAmplitudesFromCoulombIntegrals::getResiduum(
  const int iterationStep, const PTR(const FockVector<complex>) &amplitudes
) {
  throw EXCEPTION("Complex version not implemented");
}

PTR(FockVector<double>) UccsdAmplitudesFromCoulombIntegrals::getResiduum(
  const int iterationStep, const PTR(const FockVector<double>) &amplitudes
) {

  // Equations from: (It assumes hartree-fock orbitals)
  // --------------
  //John F. Stanton, Jürgen Gauss, John D. Watts, Rodney J. Bartlett.
  //A direct product decomposition approach for symmetry exploitation in
  //many‐body methods. I. Energy calculations.
  //The Journal of Chemical Physics  1991 10.1063/1.460620

  auto epsi(getTensorArgument<double>("HoleEigenEnergies"));
  auto epsa(getTensorArgument<double>("ParticleEigenEnergies"));

  // Get couloumb integrals
  auto Vijkl(getTensorArgument<double>("HHHHCoulombIntegrals"));
  auto Vabcd(getTensorArgument<double>("PPPPCoulombIntegrals"));
  auto Vabij(getTensorArgument<double>("PPHHCoulombIntegrals"));
  auto Vijka(getTensorArgument<double>("HHHPCoulombIntegrals"));
  auto Vijab(getTensorArgument<double>("HHPPCoulombIntegrals"));
  auto Viajk(getTensorArgument<double>("HPHHCoulombIntegrals"));
  auto Viajb(getTensorArgument<double>("HPHPCoulombIntegrals"));
  auto Viabc(getTensorArgument<double>("HPPPCoulombIntegrals"));
  auto Vabic(getTensorArgument<double>("PPHPCoulombIntegrals"));
  auto Viabj(getTensorArgument<double>("HPPHCoulombIntegrals"));
  auto Vaibc(getTensorArgument<double>("PHPPCoulombIntegrals"));
  auto Vijak(getTensorArgument<double>("HHPHCoulombIntegrals"));
  auto Vabci(getTensorArgument<double>("PPPHCoulombIntegrals"));

  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int vv[] = {Nv, Nv};
  int oo[] = {No, No};
  int syms[] = {NS, NS};
  CTF::Tensor<> *fab(
    new CTF::Tensor<>(2, vv, syms, *Cc4s::world, "fab")
  );
  CTF::Tensor<> *fij(
    new CTF::Tensor<>(2, oo, syms, *Cc4s::world, "fij")
  );
  CTF::Tensor<> *fia;

  if (
    isArgumentGiven("HPFockMatrix") &&
    isArgumentGiven("HHFockMatrix") &&
    isArgumentGiven("PPFockMatrix")
  ) {
    LOG(0, getAbbreviation()) << "Using non-canonical orbitals" << std::endl;
    fia = getTensorArgument<double, CTF::Tensor<> >("HPFockMatrix");
    fab = getTensorArgument<double, CTF::Tensor<> >("PPFockMatrix");
    fij = getTensorArgument<double, CTF::Tensor<> >("HHFockMatrix");
  } else {
    fia = NULL;
    (*fab)["aa"] = (*epsa)["a"];
    (*fij)["ii"] = (*epsi)["i"];
  }



  // Create T and R and intermediates
  // Read the amplitudes Tai and Tabij
  //amplitudes->get(0)
  auto Tai(amplitudes->get(0));
  Tai->set_name("Tai");
  auto Tabij(amplitudes->get(1));
  Tabij->set_name("Tabij");


  auto residuum(NEW(FockVector<double>, *amplitudes));
  *residuum *= 0.0;
  // Allocate Tensors for T2 amplitudes
  auto Rai(residuum->get(0));
  Rai->set_name("Rai");
  auto Rabij(residuum->get(1));
  Rabij->set_name("Rabij");


  // Define intermediates
  auto Fae(
    NEW(CTF::Tensor<>, 2, vv, syms, *Cc4s::world, "Fae")
  );
  auto Fmi(
    NEW(CTF::Tensor<>, 2, oo, syms, *Cc4s::world, "Fmi")
  );
  int ov[] = {No, Nv};
  auto Fme(
    NEW(CTF::Tensor<>, 2, ov, syms, *Cc4s::world, "Fme")
  );



  // Initialize intermediates

  // Equation (10)
  // TODO: Use only one Tau, since TildeTau_abij is only uysed to form Fpq
  auto Tau_abij(NEW(CTF::Tensor<>, *Tabij));
  (*Tau_abij)["abij"] += (*Tai)["ai"] * (*Tai)["bj"];
  (*Tau_abij)["abij"] += ( - 1.0 ) * (*Tai)["bi"] * (*Tai)["aj"];

  // Equation (9)
  auto TildeTau_abij(NEW(CTF::Tensor<>, *Tabij));
  (*TildeTau_abij)["abij"] += ( 0.5 ) * (*Tai)["ai"] * (*Tai)["bj"];
  (*TildeTau_abij)["abij"] += ( - 0.5 ) * (*Tai)["bi"] * (*Tai)["aj"];

  // Equation (3)
  (*Fae)["ae"]  = (*Tai)["fm"] * (*Viabc)["mafe"];
  (*Fae)["ae"] += ( - 0.5 ) * (*TildeTau_abij)["afmn"] * (*Vijab)["mnef"];

  // Equation (4)
  (*Fmi)["mi"]  = (*Tai)["en"] * (*Vijka)["mnie"];
  (*Fmi)["mi"] += (0.5) * (*TildeTau_abij)["efin"] * (*Vijab)["mnef"];

  // Equation (5) Stanton et al.
  (*Fme)["me"] = (*Tai)["fn"] * (*Vijab)["mnef"];
  if (fia) {
    (*Fme)["me"] += (*fia)["me"];
  }

  // Equation (6)
  auto Wijkl(NEW(CTF::Tensor<>, *Vijkl));
  (*Wijkl)["mnij"] += (+ 1.0) * (*Tai)["ej"] * (*Vijka)["mnie"];
  // Pij
  (*Wijkl)["mnij"] += (- 1.0) * (*Tai)["ei"] * (*Vijka)["mnje"];
  (*Wijkl)["mnij"] += (0.25) * (*Tau_abij)["efij"] * (*Vijab)["mnef"];

  // Equation (7)
  auto Wabcd(NEW(CTF::Tensor<>, *Vabcd));
  (*Wabcd)["abef"] += (- 1.0) * (*Tai)["bm"] * (*Vaibc)["amef"];
  // Pab
  (*Wabcd)["abef"] += (+ 1.0) * (*Tai)["am"] * (*Vaibc)["bmef"];
  (*Wabcd)["abef"] += (0.25) * (*Tau_abij)["abmn"] * (*Vijab)["mnef"];

  // Equation (8)
  auto Wiabj(NEW(CTF::Tensor<>, *Viabj));
  (*Wiabj)["mbej"] += (+ 1.0) * (*Tai)["fj"] * (*Viabc)["mbef"];
  (*Wiabj)["mbej"] += (- 1.0) * (*Tai)["bn"] * (*Vijak)["mnej"];
  (*Wiabj)["mbej"] += ( - 0.5 ) * (*Tabij)["fbjn"] * (*Vijab)["mnef"];
  (*Wiabj)["mbej"] +=
    ( - 1.0 ) * (*Tai)["fj"] * (*Tai)["bn"] * (*Vijab)["mnef"];

  // T1 equations:
  (*Rai)["ai"] = (*Tai)["ei"] * (*Fae)["ae"];
  if (fia) {
     (*Rai)["ai"] += (*fia)["ia"] ;
  }


  (*Rai)["ai"] += (- 1.0) * (*Tai)["am"] * (*Fmi)["mi"];
  (*Rai)["ai"] += (*Tabij)["aeim"] * (*Fme)["me"];
  (*Rai)["ai"] += (- 1.0) * (*Tai)["fn"] * (*Viajb)["naif"];
  (*Rai)["ai"] += (- 0.5) * (*Tabij)["efim"] * (*Viabc)["maef"];
  (*Rai)["ai"] += (- 0.5) * (*Tabij)["aemn"] * (*Vijak)["nmei"];

  // T2 equations:
  if (iterationStep == 0){
    LOG(1, getAbbreviation()) << "Set initial Rabij amplitudes to Vijab"
                              << std::endl;
    (*Rabij)["abij"] = (*Vijab)["ijab"];
  } else {
    (*Rabij)["abij"]  = (*Vijab)["ijab"];

    // P(ab) * Taeij ( Fbe - 0.5 Tbm Fme)
    (*Rabij)["abij"] += (1.0) * (*Tabij)["aeij"] * (*Fae)["be"];
    (*Rabij)["abij"] += (- 1.0) * (*Tabij)["beij"] * (*Fae)["ae"];
    (*Rabij)["abij"] +=
      (- 0.5) * (*Tabij)["aeij"] * (*Tai)["bm"] * (*Fme)["me"];
    (*Rabij)["abij"] +=
      (+ 0.5) * (*Tabij)["beij"] * (*Tai)["am"] * (*Fme)["me"];

    // P(ij) * Tabim ( Fmj + 0.5 Tej Fme)
    (*Rabij)["abij"] += (- 1.0) * (*Tabij)["abim"] * (*Fmi)["mj"];
    (*Rabij)["abij"] += (+ 1.0) * (*Tabij)["abjm"] * (*Fmi)["mi"];
    (*Rabij)["abij"] +=
      (- 0.5) * (*Tabij)["abim"] * (*Tai)["ej"] * (*Fme)["me"];
    (*Rabij)["abij"] +=
      (+ 0.5) * (*Tabij)["abjm"] * (*Tai)["ei"] * (*Fme)["me"];

    (*Rabij)["abij"] += (0.5) * (*Tau_abij)["abmn"] * (*Wijkl)["mnij"];
    (*Rabij)["abij"] += (0.5) * (*Tau_abij)["efij"] * (*Wabcd)["abef"];

    // P-ij * P-ab
    (*Rabij)["abij"] += (  1.0) * (*Tabij)["aeim"] * (*Wiabj)["mbej"];
    // -Pij
    (*Rabij)["abij"] += (- 1.0) * (*Tabij)["aejm"] * (*Wiabj)["mbei"];
    // -Pab
    (*Rabij)["abij"] += (- 1.0) * (*Tabij)["beim"] * (*Wiabj)["maej"];
    //  Pij * Pab
    (*Rabij)["abij"] += (  1.0) * (*Tabij)["bejm"] * (*Wiabj)["maei"];

    // P-ij * P-ab
    (*Rabij)["abij"] +=
      (- 1.0) * (*Tai)["ei"] * (*Tai)["am"] * (*Viabj)["mbej"];
    // +Pij
    (*Rabij)["abij"] +=
      (+ 1.0) * (*Tai)["ej"] * (*Tai)["am"] * (*Viabj)["mbei"];
    // +Pab
    (*Rabij)["abij"] +=
      (+ 1.0) * (*Tai)["ei"] * (*Tai)["bm"] * (*Viabj)["maej"];
    //  - Pij * Pab
    (*Rabij)["abij"] +=
      (- 1.0) * (*Tai)["ej"] * (*Tai)["bm"] * (*Viabj)["maei"];

    (*Rabij)["abij"] += (  1.0) * (*Tai)["ei"] * (*Vabci)["abej"];
    // - Pij
    (*Rabij)["abij"] += (- 1.0) * (*Tai)["ej"] * (*Vabci)["abei"];

    (*Rabij)["abij"] += (- 1.0) * (*Tai)["am"] * (*Viajk)["mbij"];
    // + Pab
    (*Rabij)["abij"] += (+ 1.0) * (*Tai)["bm"] * (*Viajk)["maij"];

  }

  return residuum;
}



