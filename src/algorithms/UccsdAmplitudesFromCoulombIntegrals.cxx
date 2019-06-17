#include <algorithms/UccsdAmplitudesFromCoulombIntegrals.hpp>
#include <unistd.h>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <math/MathFunctions.hpp>
#include <math/ComplexTensor.hpp>
#include <math/RandomTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <util/RangeParser.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>

using namespace cc4s;
using namespace tcc;

ALGORITHM_REGISTRAR_DEFINITION(UccsdAmplitudesFromCoulombIntegrals);

UccsdAmplitudesFromCoulombIntegrals::UccsdAmplitudesFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): ClusterSinglesDoublesAlgorithm(argumentList) {}


UccsdAmplitudesFromCoulombIntegrals::~UccsdAmplitudesFromCoulombIntegrals() {
}

void UccsdAmplitudesFromCoulombIntegrals::run() {
  usingIntermediates = (bool) getIntegerArgument("intermediates", 1);
  if (! usingIntermediates ) {
    LOG(0, getAbbreviation()) <<
      "Not using intermediates, the code will be much slower."
    << std::endl;
  }
  ClusterSinglesDoublesAlgorithm::run();
}


PTR(FockVector<cc4s::complex>) UccsdAmplitudesFromCoulombIntegrals::getResiduum(
  const int iterationStep, const PTR(const FockVector<complex>) &amplitudes
) {
  if (iterationStep == 0){
    LOG(1, getAbbreviation()) <<
       "WARNING: Using complex version of Uccsd" << std::endl;
    LOG(1, getAbbreviation()) <<
       "WARNING: Complex version is not tested." << std::endl;
  }
  return getResiduumTemplate<complex>(iterationStep, amplitudes);
}

PTR(FockVector<double>) UccsdAmplitudesFromCoulombIntegrals::getResiduum(
  const int iterationStep, const PTR(const FockVector<double>) &amplitudes
) {
  return getResiduumTemplate<double>(iterationStep, amplitudes);
}

template <typename F>
PTR(FockVector<F>) UccsdAmplitudesFromCoulombIntegrals::getResiduumTemplate(
  const int iterationStep, const PTR(const FockVector<F>) &amplitudes
) {
  // Equations from:
  // --------------
  //John F. Stanton, Jürgen Gauss, John D. Watts, Rodney J. Bartlett.
  //A direct product decomposition approach for symmetry exploitation in
  //many‐body methods. I. Energy calculations.
  //The Journal of Chemical Physics  1991 doi:10.1063/1.460620

  CTF::Tensor<double> *epsi(
    getTensorArgument<double, CTF::Tensor<double> >("HoleEigenEnergies")
  );

  CTF::Tensor<double> *epsa(
    getTensorArgument<double, CTF::Tensor<double> >("ParticleEigenEnergies")
  );

  // Get couloumb integrals
  auto Vijkl(getTensorArgument<F, CTF::Tensor<F> >("HHHHCoulombIntegrals"));
  auto Vabcd(getTensorArgument<F, CTF::Tensor<F> >("PPPPCoulombIntegrals"));
  auto Vijka(getTensorArgument<F, CTF::Tensor<F> >("HHHPCoulombIntegrals"));
  auto Vijab(getTensorArgument<F, CTF::Tensor<F> >("HHPPCoulombIntegrals"));
  auto Viajk(getTensorArgument<F, CTF::Tensor<F> >("HPHHCoulombIntegrals"));
  auto Viajb(getTensorArgument<F, CTF::Tensor<F> >("HPHPCoulombIntegrals"));
  auto Viabc(getTensorArgument<F, CTF::Tensor<F> >("HPPPCoulombIntegrals"));
  auto Vabij(getTensorArgument<F, CTF::Tensor<F> >("PPHHCoulombIntegrals"));
  auto Vabic(getTensorArgument<F, CTF::Tensor<F> >("PPHPCoulombIntegrals"));
  auto Viabj(getTensorArgument<F, CTF::Tensor<F> >("HPPHCoulombIntegrals"));
  auto Vaibc(getTensorArgument<F, CTF::Tensor<F> >("PHPPCoulombIntegrals"));
  auto Vijak(getTensorArgument<F, CTF::Tensor<F> >("HHPHCoulombIntegrals"));
  auto Vabci(getTensorArgument<F, CTF::Tensor<F> >("PPPHCoulombIntegrals"));

  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int vv[] = {Nv, Nv};
  int oo[] = {No, No};
  int syms[] = {NS, NS};
  CTF::Tensor<F> *fab(
    new CTF::Tensor<F>(2, vv, syms, *Cc4s::world, "fab")
  );
  CTF::Tensor<F> *fij(
    new CTF::Tensor<F>(2, oo, syms, *Cc4s::world, "fij")
  );
  CTF::Tensor<F> *fia;

  if (
    isArgumentGiven("HPFockMatrix") &&
    isArgumentGiven("HHFockMatrix") &&
    isArgumentGiven("PPFockMatrix")
  ) {
    if (iterationStep == 0){
    LOG(0, getAbbreviation()) << "Using non-canonical orbitals" << std::endl;
    }
    fia = getTensorArgument<F, CTF::Tensor<F> >("HPFockMatrix");
    fab = getTensorArgument<F, CTF::Tensor<F> >("PPFockMatrix");
    fij = getTensorArgument<F, CTF::Tensor<F> >("HHFockMatrix");
  } else {
    fia = NULL;
    CTF::Transform<double, F>(
      std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }
      )
    ) (
      (*epsi)["i"], (*fij)["ii"]
    );
    CTF::Transform<double, F>(
      std::function<void(double, F &)>(
        [](double eps, F &f) { f = eps; }
      )
    ) (
      (*epsa)["a"], (*fab)["aa"]
    );
  }

  // Create T and R and intermediates
  // Read the amplitudes Tai and Tabij
  //amplitudes->get(0)
  auto Tai(amplitudes->get(0));
  Tai->set_name("Tai");
  auto Tabij(amplitudes->get(1));
  Tabij->set_name("Tabij");


  auto residuum(NEW(FockVector<F>, *amplitudes));
  *residuum *= 0.0;
  // Allocate Tensors for T2 amplitudes
  auto Rai(residuum->get(0));
  Rai->set_name("Rai");
  auto Rabij(residuum->get(1));
  Rabij->set_name("Rabij");

  if (usingIntermediates) {

    // Define intermediates
    auto Fae(
      NEW(CTF::Tensor<F>, 2, vv, syms, *Cc4s::world, "Fae")
    );
    auto Fmi(
      NEW(CTF::Tensor<F>, 2, oo, syms, *Cc4s::world, "Fmi")
    );
    int ov[] = {No, Nv};
    auto Fme(
      NEW(CTF::Tensor<F>, 2, ov, syms, *Cc4s::world, "Fme")
    );

    // Initialize intermediates

    // Equation (10)
    auto Tau_abij(NEW(CTF::Tensor<F>, *Tabij));
    (*Tau_abij)["abij"] += (*Tai)["ai"] * (*Tai)["bj"];
    (*Tau_abij)["abij"] += ( - 1.0 ) * (*Tai)["bi"] * (*Tai)["aj"];

    // Equation (9)
    auto TildeTau_abij(NEW(CTF::Tensor<F>, *Tabij));
    (*TildeTau_abij)["abij"] += ( 0.5 ) * (*Tai)["ai"] * (*Tai)["bj"];
    (*TildeTau_abij)["abij"] += ( - 0.5 ) * (*Tai)["bi"] * (*Tai)["aj"];

    // Equation (3)
    (*Fae)["ae"]  = (*fab)["ae"];
    (*Fae)["aa"] += (-1.0) * (*fab)["aa"];
    if (fia) {
      (*Fae)["aa"] += (-0.5) * (*fia)["me"] * (*Tai)["am"];
    }
    (*Fae)["ae"] += (*Tai)["fm"] * (*Viabc)["mafe"];
    (*Fae)["ae"] += ( - 0.5 ) * (*TildeTau_abij)["afmn"] * (*Vijab)["mnef"];

    // Equation (4)
    (*Fmi)["mi"]  = (*fij)["mi"];
    (*Fmi)["ii"] += (-1.0) * (*fij)["ii"];
    if (fia) {
      (*Fmi)["mi"] += (+0.5) * (*fia)["me"] * (*Tai)["ei"];
    }
    (*Fmi)["mi"] += (*Tai)["en"] * (*Vijka)["mnie"];
    (*Fmi)["mi"] += (0.5) * (*TildeTau_abij)["efin"] * (*Vijab)["mnef"];

    // Equation (5) Stanton et al.
    (*Fme)["me"] = (*Tai)["fn"] * (*Vijab)["mnef"];
    if (fia) {
      (*Fme)["me"] += (*fia)["me"];
    }

    // Equation (6)
    auto Wijkl(NEW(CTF::Tensor<F>, *Vijkl));
    (*Wijkl)["mnij"] += (+ 1.0) * (*Tai)["ej"] * (*Vijka)["mnie"];
    // Pij
    (*Wijkl)["mnij"] += (- 1.0) * (*Tai)["ei"] * (*Vijka)["mnje"];
    (*Wijkl)["mnij"] += (0.25) * (*Tau_abij)["efij"] * (*Vijab)["mnef"];

    // Equation (7)
    auto Wabcd(NEW(CTF::Tensor<F>, *Vabcd));
    (*Wabcd)["abef"] += (- 1.0) * (*Tai)["bm"] * (*Vaibc)["amef"];
    // Pab
    (*Wabcd)["abef"] += (+ 1.0) * (*Tai)["am"] * (*Vaibc)["bmef"];
    (*Wabcd)["abef"] += (0.25) * (*Tau_abij)["abmn"] * (*Vijab)["mnef"];

    // Equation (8)
    auto Wiabj(NEW(CTF::Tensor<F>, *Viabj));
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

      return residuum;

    }

    //fbd136af1cf7b395e33a6a2040e36b92d7a18e27  -
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

    //dda7881e81095a6a21011833e88cc972ff890456  -
    //102d94b283918c672a8c9ce73ef689f2f8740661  -
    (*Rabij)["abij"] += (0.5) * (*Tau_abij)["efij"] * (*Wabcd)["abef"];

    //ec0d590a53887b24c495ca90df46ee5782c62515  -
    //e7b0e03e493d603550574c6f4015d16f994eb184  -
    //cff259570ee87e7b824936881b03dc0dc3e55c80  -
    //bcaaa0ebd5951466bf1b3592a990c241c26c7b4e  -
    //aa552987b73a1b0556efa1ab1bb7edf26d5ed58d  -
    // P-ij * P-ab
    (*Rabij)["abij"] += (  1.0) * (*Tabij)["aeim"] * (*Wiabj)["mbej"];
    // -Pij
    (*Rabij)["abij"] += (- 1.0) * (*Tabij)["aejm"] * (*Wiabj)["mbei"];
    // -Pab
    (*Rabij)["abij"] += (- 1.0) * (*Tabij)["beim"] * (*Wiabj)["maej"];
    //  Pij * Pab
    (*Rabij)["abij"] += (  1.0) * (*Tabij)["bejm"] * (*Wiabj)["maei"];

    //896383ad3db1c77dc734c836141429c67784e118  -
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

    //e035d48a19d337a004e70549c1a490cba713f419  -
    (*Rabij)["abij"] += (  1.0) * (*Tai)["ei"] * (*Vabci)["abej"];
    // - Pij
    (*Rabij)["abij"] += (- 1.0) * (*Tai)["ej"] * (*Vabci)["abei"];

    //b7a00eb9acb88bcf1e3fc73992e621beb09b39f2  -
    (*Rabij)["abij"] += (- 1.0) * (*Tai)["am"] * (*Viajk)["mbij"];
    // + Pab
    (*Rabij)["abij"] += (+ 1.0) * (*Tai)["bm"] * (*Viajk)["maij"];

    // End intermediates

  } else {

    (*Rabij)["abij"]  = 0;
    (*Rai)["ai"]  = 0;

    // Equations from hirata
    //
    // Singles <bi| e^-T H e^T|0> = 0
    if (fia) {
      //(*Rai)["bi"] += ( + 1.0  ) * (*fai)["bi"];
      (*Rai)["bi"] += ( + 1.0  ) * (*fia)["ib"];
      (*Rai)["bi"] += ( + 1.0  ) * (*fia)["kd"] * (*Tabij)["dbki"];
      (*Rai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["bl"] * (*fia)["lc"];
    }

    //These are the residum equations
    //(*Rai)["bi"] += ( + 1.0  ) * (*fab)["bc"] * (*Tai)["ci"];
    //(*Rai)["bi"] += ( - 1.0  ) * (*fij)["ki"] * (*Tai)["bk"];

    (*Rai)["bi"] += ( - 1.0  ) * (*Tai)["cl"] * (*Viajb)["lbic"];
    (*Rai)["bi"] += ( + 0.5  ) * (*Tabij)["cblm"] * (*Vijka)["lmic"];
    (*Rai)["bi"] += ( + 0.5  ) * (*Tabij)["cdmi"] * (*Viabc)["mbcd"];
    (*Rai)["bi"] += ( - 1.0  ) * (*Tai)["bk"] * (*Tai)["dm"] * (*Vijka)["kmid"];
    (*Rai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["dm"] * (*Viabc)["mbcd"];
    (*Rai)["bi"] += ( - 0.5  ) * (*Tai)["fi"] * (*Tabij)["cblm"] * (*Vijab)["lmcf"];
    (*Rai)["bi"] += ( - 0.5  ) * (*Tai)["bn"] * (*Tabij)["cdmi"] * (*Vijab)["mncd"];
    (*Rai)["bi"] += ( + 1.0  ) * (*Tabij)["cbli"] * (*Tai)["en"] * (*Vijab)["lnce"];
    (*Rai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["bl"] * (*Tai)["en"] * (*Vijab)["lnce"];

    // Doubles <cdij| e^-T H e^T|0> = 0
    //fbd136af1cf7b395e33a6a2040e36b92d7a18e27  -
    (*Rabij)["cdij"]  = ( + 1.0  ) * (*Vabij)["cdij"];

    if (fia) {
      (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["cdmj"] * (*fia)["mf"] * (*Tai)["fi"];
      (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["cdmi"] * (*fia)["mf"] * (*Tai)["fj"];
      (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["fcij"] * (*fia)["mf"] * (*Tai)["dm"];
      (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["fdij"] * (*fia)["mf"] * (*Tai)["cm"];
    }

    //These are the residum equations
    //(*Rabij)["cdij"] += ( - 1.0  ) * (*fij)["mi"] * (*Tabij)["cdmj"];
    //(*Rabij)["cdij"] += ( + 1.0  ) * (*fij)["mj"] * (*Tabij)["cdmi"];
    //(*Rabij)["cdij"] += ( - 1.0  ) * (*fab)["de"] * (*Tabij)["ecij"];
    //(*Rabij)["cdij"] += ( + 1.0  ) * (*fab)["ce"] * (*Tabij)["edij"];

    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Vijkl)["mnij"];

    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijkl)["mnij"];

    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gj"] * (*Vijka)["mnig"];
    (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gi"] * (*Vijka)["mnjg"];

    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["cdmj"] * (*Tai)["fo"] * (*Vijka)["moif"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["cdmi"] * (*Tai)["fo"] * (*Vijka)["mojf"];

    (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabij)["efij"] * (*Tai)["co"] * (*Viabc)["odef"];
    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["efij"] * (*Tai)["do"] * (*Viabc)["ocef"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecij"] * (*Tai)["fo"] * (*Viabc)["odef"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edij"] * (*Tai)["fo"] * (*Viabc)["ocef"];

    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["edij"] * (*Tabij)["fcop"] * (*Vijab)["opef"];
    (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabij)["ecij"] * (*Tabij)["fdop"] * (*Vijab)["opef"];

    (*Rabij)["cdij"] += ( + 0.25  ) * (*Tabij)["efij"] * (*Tabij)["cdop"] * (*Vijab)["opef"];

    (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabij)["cdmi"] * (*Tabij)["fgpj"] * (*Vijab)["mpfg"];
    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["cdmj"] * (*Tabij)["fgpi"] * (*Vijab)["mpfg"];

    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijka)["noie"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijka)["noje"];

    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Viabc)["odef"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["do"] * (*Viabc)["ocef"];

    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tabij)["cdop"] * (*Vijab)["opef"];

    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["cdpj"] * (*Tai)["ei"] * (*Tai)["fo"] * (*Vijab)["opef"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["cdpi"] * (*Tai)["ej"] * (*Tai)["fo"] * (*Vijab)["opef"];

    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Tabij)["ghij"] * (*Vijab)["mngh"];

    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["hcij"] * (*Tai)["dm"] * (*Tai)["fo"] * (*Vijab)["mofh"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["hdij"] * (*Tai)["cm"] * (*Tai)["fo"] * (*Vijab)["mofh"];

    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Tai)["dp"] * (*Vijab)["opef"];

    //dda7881e81095a6a21011833e88cc972ff890456  -
    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["efij"] * (*Vabcd)["cdef"];
    //102d94b283918c672a8c9ce73ef689f2f8740661  -
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Vabcd)["cdef"];

    //ec0d590a53887b24c495ca90df46ee5782c62515  -
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ecni"] * (*Viajb)["ndje"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecnj"] * (*Viajb)["ndie"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["edni"] * (*Viajb)["ncje"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ednj"] * (*Viajb)["ncie"];
    //e7b0e03e493d603550574c6f4015d16f994eb184  -
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Tai)["gj"] * (*Viabc)["ndeg"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["gi"] * (*Viabc)["ndeg"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Tai)["gj"] * (*Viabc)["nceg"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ednj"] * (*Tai)["gi"] * (*Viabc)["nceg"];
    //cff259570ee87e7b824936881b03dc0dc3e55c80  -
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Tai)["do"] * (*Vijka)["noje"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["do"] * (*Vijka)["noie"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Tai)["co"] * (*Vijka)["noje"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ednj"] * (*Tai)["co"] * (*Vijka)["noie"];
    //bcaaa0ebd5951466bf1b3592a990c241c26c7b4e  -
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Tabij)["gcpj"] * (*Vijab)["npeg"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Tabij)["gdpj"] * (*Vijab)["npeg"];
    //aa552987b73a1b0556efa1ab1bb7edf26d5ed58d  -
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["gcpi"] * (*Tai)["ej"] * (*Tai)["dn"] * (*Vijab)["npeg"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["gcpj"] * (*Tai)["ei"] * (*Tai)["dn"] * (*Vijab)["npeg"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["gdpi"] * (*Tai)["ej"] * (*Tai)["cn"] * (*Vijab)["npeg"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["gdpj"] * (*Tai)["ei"] * (*Tai)["cn"] * (*Vijab)["npeg"];

    //896383ad3db1c77dc734c836141429c67784e118  -
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Viajb)["ndje"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Viajb)["ndie"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Viajb)["ncje"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Viajb)["ncie"];

    //e035d48a19d337a004e70549c1a490cba713f419  -
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Vabic)["cdie"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Vabic)["cdje"];

    //b7a00eb9acb88bcf1e3fc73992e621beb09b39f2  -
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["cm"] * (*Viajk)["mdij"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["dm"] * (*Viajk)["mcij"];

  }

  return residuum;
}



