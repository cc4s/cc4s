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

double UccsdAmplitudesFromCoulombIntegrals::calculateEnergy() {
  // allocate energy
  CTF::Scalar<double> energy(*Vabij->wrld);
  energy.set_name("energy");

  // get amplitudes from the mixers
  CTF::Tensor<> *Tai(&TaiMixer->getNext());
  Tai->set_name("Tai");
  CTF::Tensor<> *Tabij(&TabijMixer->getNext());
  Tabij->set_name("Tabij");

  energy[""] += ( + 0.25  ) * (*Tabij)["abkl"] * (*Vijab)["klab"];
  energy[""] += ( + 0.5  ) * (*Tai)["aj"] * (*Tai)["cl"] * (*Vijab)["jlac"];

  double e(energy.get_val());
  LOG(0, "UCCSD") << "e=" << e << std::endl;

  return e;
}


void UccsdAmplitudesFromCoulombIntegrals::run() {
  LOG(0, "UCCSD") << "Initializing data" << std::endl;
  epsi = getTensorArgument<double, CTF::Tensor<> >("HoleEigenEnergies");
  epsa = getTensorArgument<double, CTF::Tensor<> >("ParticleEigenEnergies");

  int Nv(epsa->lens[0]), No(epsi->lens[0]);

  // Get couloumb integrals (these shoul not be antisymmetrized)
  Vijkl = getTensorArgument<double, CTF::Tensor<> >("HHHHCoulombIntegrals");
  Vabcd = getTensorArgument<double, CTF::Tensor<> >("PPPPCoulombIntegrals");
  Vabij = getTensorArgument<double, CTF::Tensor<> >("PPHHCoulombIntegrals");
  Vijka = getTensorArgument<double, CTF::Tensor<> >("HHHPCoulombIntegrals");
  Vaibj = getTensorArgument<double, CTF::Tensor<> >("PHPHCoulombIntegrals");
  Vabci = getTensorArgument<double, CTF::Tensor<> >("PPPHCoulombIntegrals");

  int syms[] = {NS, NS, NS, NS};

  //  Vijab
  int oovv[] = { No, No, Nv, Nv };
  Vijab =  new CTF::Tensor<>(4, oovv, syms, *Cc4s::world, "Vijab");
  (*Vijab)["ijab"] =  (*Vabij)["abij"] - (*Vabij)["abji"];
  conjugate(*Vijab);

  //  Viajk
  int ovoo[] = { No, Nv, No, No };
  Viajk =  new CTF::Tensor<>(4, ovoo, syms, *Cc4s::world, "Viajk");
  (*Viajk)["iajk"] =  (*Vijka)["jkia"]  - (*Vijka)["kjia"];
  conjugate(*Viajk);

  // Viajb
  int ovov[] = { No, Nv, No, Nv };
  int voov[] = { Nv, No, No, Nv };
  CTF::Tensor<> Vaijb(4, voov, syms, *Cc4s::world, "Vaijb");
  // Assumes real orbitals
  Vaijb["aijb"] = (*Vabij)["abji"];
  Viajb =  new CTF::Tensor<>(4, ovov, syms, *Cc4s::world, "Viajb");
  (*Viajb)["iajb"] = ( - 1.0 ) * (Vaijb)["aijb"];
  (*Viajb)["iajb"] +=  (*Vaibj)["aibj"];

  // Viabc
  int ovvv[] = { No, Nv, Nv, Nv };
  Viabc =  new CTF::Tensor<>(4, ovvv, syms, *Cc4s::world, "Viabc");
  (*Viabc)["iabc"] =  (*Vabci)["abci"];
  conjugate(*Viabc);
  (*Viabc)["iabc"] -= (*Vabci)["acbi"];

  // Vabic
  int vvov[] = { Nv, Nv, No, Nv };
  Vabic =  new CTF::Tensor<>(4, vvov, syms, *Cc4s::world, "Vabic");
  (*Vabic)["abic"] =  (*Vabci)["baci"] - (*Vabci)["abci"];

  // Antisymmetrize integrals that are read in
  (*Vijkl)["ijkl"] -= (*Vijkl)["ijlk"];
  (*Vabcd)["abcd"] -= (*Vabcd)["abdc"];
  (*Vijka)["ijka"] -= (*Vijka)["jika"];
  (*Vaibj)["aibj"] -= (*Vabij)["baij"];
  (*Vabci)["abci"] -= (*Vabci)["baci"];
  (*Vabij)["abij"] -= (*Vabij)["abji"];

  ClusterSinglesDoublesAlgorithm::run();
}
void UccsdAmplitudesFromCoulombIntegrals::iterate(int iterationStep) {

  int Nv(epsa->lens[0]), No(epsi->lens[0]);

  // Create T and R and intermediates
  // Read the amplitudes Tai and Tabij
  CTF::Tensor<> *Tai(&TaiMixer->getNext());
  Tai->set_name("Tai");
  CTF::Tensor<> *Tabij(&TabijMixer->getNext());
  Tabij->set_name("Tabij");

  // Allocate Tensors for T2 amplitudes
  CTF::Tensor<> Rai(*Tai);
  Rai.set_name("Rai");
  CTF::Tensor<> Rabij(*Vabij);
  Rabij.set_name("Rabij");
  //Tai->print();
  //Tabij->print();

  // kinetic terms
  int oneBodySyms[] = {NS, NS};
  int vv[] = {Nv, Nv};
  CTF::Tensor<> *Fab( new CTF::Tensor<>(2, vv, oneBodySyms, *Cc4s::world, "Fab") );
  int oo[] = {No, No};
  CTF::Tensor<> *Fij( new CTF::Tensor<>(2, oo, oneBodySyms, *Cc4s::world, "Fij") );

  (*Fab)["aa"] = (*epsa)["a"];
  (*Fij)["ii"] = (*epsi)["i"];

  Rai["bi"]  = 0.0;
  Rai["bi"] += ( - 1.0  ) * (*Fij)["ki"] * (*Tai)["bk"];
  Rai["bi"] += ( + 1.0  ) * (*Fab)["bc"] * (*Tai)["ci"];
  Rai["bi"] += ( - 1.0  ) * (*Tai)["cl"] * (*Viajb)["lbic"];
  Rai["bi"] += ( + 0.5  ) * (*Tabij)["cblm"] * (*Vijka)["lmic"];
  Rai["bi"] += ( + 0.5  ) * (*Tabij)["cdmi"] * (*Viabc)["mbcd"];
  Rai["bi"] += ( - 1.0  ) * (*Tai)["bk"] * (*Tai)["dm"] * (*Vijka)["kmid"];
  Rai["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["dm"] * (*Viabc)["mbcd"];
  Rai["bi"] += ( - 0.5  ) * (*Tabij)["cblm"] * (*Tai)["fi"] * (*Vijab)["lmcf"];
  Rai["bi"] += ( - 0.5  ) * (*Tabij)["cdmi"] * (*Tai)["bn"] * (*Vijab)["mncd"];
  Rai["bi"] += ( + 1.0  ) * (*Tabij)["cbli"] * (*Tai)["en"] * (*Vijab)["lnce"];
  Rai["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["bl"] * (*Tai)["en"] * (*Vijab)["lnce"];

  singlesAmplitudesFromResiduum(Rai);
  TaiMixer->append(Rai);

  if (iterationStep == 0){
    LOG(1, getAbbreviation()) << "Set initial Rabij amplitudes to Vabij"
                              << std::endl;
    Rabij["abij"] = (*Vabij)["abij"];
  } else {
    Rabij["cdij"]  = 0.0;
    Rabij["cdij"] += ( + 1.0  ) * (*Vabij)["cdij"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tai)["cm"] * (*Viajk)["mdij"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tai)["dm"] * (*Viajk)["mcij"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Vabic)["cdie"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Vabic)["cdje"];
    Rabij["cdij"] += ( - 1.0  ) * (*Fij)["mi"] * (*Tabij)["cdmj"];
    Rabij["cdij"] += ( + 1.0  ) * (*Fij)["mj"] * (*Tabij)["cdmi"];
    Rabij["cdij"] += ( - 1.0  ) * (*Fab)["de"] * (*Tabij)["ecij"];
    Rabij["cdij"] += ( + 1.0  ) * (*Fab)["ce"] * (*Tabij)["edij"];
    Rabij["cdij"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Vijkl)["mnij"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tabij)["ecnj"] * (*Viajb)["ndie"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tabij)["ednj"] * (*Viajb)["ncie"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tabij)["ecni"] * (*Viajb)["ndje"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tabij)["edni"] * (*Viajb)["ncje"];
    Rabij["cdij"] += ( + 0.5  ) * (*Tabij)["efij"] * (*Vabcd)["cdef"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijkl)["mnij"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Viajb)["ndie"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Viajb)["ncie"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Viajb)["ndje"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Viajb)["ncje"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Vabcd)["cdef"];
    Rabij["cdij"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gj"] * (*Vijka)["mnig"];
    Rabij["cdij"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gi"] * (*Vijka)["mnjg"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["do"] * (*Vijka)["noie"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tabij)["ednj"] * (*Tai)["co"] * (*Vijka)["noie"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Tai)["do"] * (*Vijka)["noje"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Tai)["co"] * (*Vijka)["noje"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tabij)["cdmj"] * (*Tai)["fo"] * (*Vijka)["moif"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tabij)["cdmi"] * (*Tai)["fo"] * (*Vijka)["mojf"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["gi"] * (*Viabc)["ndeg"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tabij)["ednj"] * (*Tai)["gi"] * (*Viabc)["nceg"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Tai)["gj"] * (*Viabc)["ndeg"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Tai)["gj"] * (*Viabc)["nceg"];
    Rabij["cdij"] += ( - 0.5  ) * (*Tabij)["efij"] * (*Tai)["co"] * (*Viabc)["odef"];
    Rabij["cdij"] += ( + 0.5  ) * (*Tabij)["efij"] * (*Tai)["do"] * (*Viabc)["ocef"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tabij)["ecij"] * (*Tai)["fo"] * (*Viabc)["odef"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tabij)["edij"] * (*Tai)["fo"] * (*Viabc)["ocef"];
    Rabij["cdij"] += ( + 0.5  ) * (*Tabij)["edij"] * (*Tabij)["fcop"] * (*Vijab)["opef"];
    Rabij["cdij"] += ( - 0.5  ) * (*Tabij)["ecij"] * (*Tabij)["fdop"] * (*Vijab)["opef"];
    Rabij["cdij"] += ( + 0.25  ) * (*Tabij)["efij"] * (*Tabij)["cdop"] * (*Vijab)["opef"];
    Rabij["cdij"] += ( - 0.5  ) * (*Tabij)["cdmi"] * (*Tabij)["fgpj"] * (*Vijab)["mpfg"];
    Rabij["cdij"] += ( + 0.5  ) * (*Tabij)["cdmj"] * (*Tabij)["fgpi"] * (*Vijab)["mpfg"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Tabij)["gcpj"] * (*Vijab)["npeg"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Tabij)["gdpj"] * (*Vijab)["npeg"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijka)["noie"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijka)["noje"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Viabc)["odef"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["do"] * (*Viabc)["ocef"];
    Rabij["cdij"] += ( + 0.5  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tabij)["cdop"] * (*Vijab)["opef"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Tabij)["gcpj"] * (*Vijab)["npeg"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Tabij)["gdpj"] * (*Vijab)["npeg"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Tabij)["gcpi"] * (*Vijab)["npeg"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Tabij)["gdpi"] * (*Vijab)["npeg"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fo"] * (*Tabij)["cdpj"] * (*Vijab)["opef"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["fo"] * (*Tabij)["cdpi"] * (*Vijab)["opef"];
    Rabij["cdij"] += ( + 0.5  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Tabij)["ghij"] * (*Vijab)["mngh"];
    Rabij["cdij"] += ( - 1.0  ) * (*Tai)["dm"] * (*Tai)["fo"] * (*Tabij)["hcij"] * (*Vijab)["mofh"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["fo"] * (*Tabij)["hdij"] * (*Vijab)["mofh"];
    Rabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Tai)["dp"] * (*Vijab)["opef"];
  }

  // Calculate the amplitudes from the residuum
  doublesAmplitudesFromResiduum(Rabij);
  // Append amplitudes to the mixer
  TabijMixer->append(Rabij);


}



