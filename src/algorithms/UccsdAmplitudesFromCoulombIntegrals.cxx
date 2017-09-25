#include <algorithms/UccsdAmplitudesFromCoulombIntegrals.hpp>
#include <unistd.h>
#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <math/MathFunctions.hpp>
#include <math/RandomTensor.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>

using namespace cc4s;
using namespace tcc;

ALGORITHM_REGISTRAR_DEFINITION(UccsdAmplitudesFromCoulombIntegrals);

// TODO: Study the requirements to treat the problem with real numbers or
// complex

UccsdAmplitudesFromCoulombIntegrals::UccsdAmplitudesFromCoulombIntegrals(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

UccsdAmplitudesFromCoulombIntegrals::~UccsdAmplitudesFromCoulombIntegrals() {
}

void UccsdAmplitudesFromCoulombIntegrals::run() {
  typedef CTF::Tensor<> T;

  // Get orbital energies
  T *epsi(getTensorArgument<double, T>("HoleEigenEnergies"));
  T *epsa(getTensorArgument<double, T>("ParticleEigenEnergies"));

  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int totalDimension(1 + Nv * No + No * No * Nv * Nv);
  LOG(1, "Uccsd_AMPS") << "Nv " << Nv << std::endl;
  LOG(1, "Uccsd_AMPS") << "No " << No << std::endl;
  LOG(1, "Uccsd_AMPS") << "Problem dimension " << totalDimension << std::endl;

  // Get couloumb integrals (these shoul not be antisymmetrized)
  T *Vijkl(getTensorArgument<double, T>("HHHHCoulombIntegrals"));
  T *Vabcd(getTensorArgument<double, T>("PPPPCoulombIntegrals")); 

  T *Vabij(getTensorArgument<double, T>("PPHHCoulombIntegrals"));
  // T *Vijab(getTensorArgument<double, T>("HHPPCoulombIntegrals")); // swap PPHH (done)

  T *Vijka(getTensorArgument<double, T>("HHHPCoulombIntegrals"));
  // T *Viajk(getTensorArgument<double, T>("HPHHCoulombIntegrals")); // swap HHHP (done)

  T *Vaibj(getTensorArgument<double, T>("PHPHCoulombIntegrals")); // not in eqs
  //T *Viajb(getTensorArgument<double, T>("HPHPCoulombIntegrals")); // swap PHPH (done)

  T *Vabci(getTensorArgument<double, T>("PPPHCoulombIntegrals")); // not in eqs
  //T *Viabc(getTensorArgument<double, T>("HPPPCoulombIntegrals")); // swap PPPH (done)
  //T *Vabic(getTensorArgument<double, T>("PPHPCoulombIntegrals")); // swap PPPH (done)

  int syms[] = {NS, NS, NS, NS};

  LOG(1, "Uccsd_AMPS") << "Antisymmetrizing Vpqrs " << std::endl;

  //  Vijab
  int oovv[] = { No, No, Nv, Nv };
  T *Vijab(
    new T(4, oovv, syms, *Cc4s::world, "Vijab")
  );
  (*Vijab)["ijab"] =  (*Vabij)["abij"] - (*Vabij)["abji"];

  //  Viajk
  int ovoo[] = { No, Nv, No, No };
  T *Viajk(
    new T(4, ovoo, syms, *Cc4s::world, "Viajk")
  );
  (*Viajk)["iajk"] =  (*Vijka)["ijka"]  - (*Vijka)["ikja"];

  // Viajb
  int ovov[] = { No, Nv, No, Nv };
  T *Viajb(
    new T(4, ovov, syms, *Cc4s::world, "Viajb")
  );
  (*Viajb)["iajb"] =  (*Vaibj)["aibj"] - (*Vabij)["abji"];

  // Viabc
  int ovvv[] = { No, Nv, Nv, Nv };
  T *Viabc(
    new T(4, ovvv, syms, *Cc4s::world, "Viabc")
  );
  (*Viabc)["iabc"] =  (*Vabci)["abci"] - (*Vabci)["acbi"];

  // Vabic
  int vvov[] = { Nv, Nv, No, Nv };
  T *Vabic(
    new T(4, vvov, syms, *Cc4s::world, "Vabic")
  );
  (*Vabic)["abic"] =  (*Vabci)["abci"] - (*Vabci)["baci"];

  // Antisymmetrize integrals that are read in
  (*Vijkl)["ijkl"] -= (*Vijkl)["ijlk"];
  (*Vabcd)["abcd"] -= (*Vabcd)["abdc"];
  (*Vijka)["ijka"] -= (*Vijka)["jika"];
  (*Vaibj)["aibj"] -= (*Vabij)["baij"];
  (*Vabci)["abci"] -= (*Vabci)["baci"];
  (*Vabij)["abij"] -= (*Vabij)["abji"];


  // Create T and R and intermediates
  int oneBodySyms[] = {NS, NS};
  int vo[] = {Nv, No};
  int vvoo[] = {Nv, Nv, No, No};
  T *Tai( new T(2, vo, oneBodySyms, *Cc4s::world, "Tai") );
  T *Tabij(new T(4, vvoo, syms, *Cc4s::world, "Tabij") );
  T *Rai( new T(2, vo, oneBodySyms, *Cc4s::world, "Rai") );
  T *Rabij(new T(4, vvoo, syms, *Cc4s::world, "Rabij") );

  // kinetic terms
  int vv[] = {Nv, Nv};
  T *Fab( new T(2, vv, oneBodySyms, *Cc4s::world, "Fab") );
  int oo[] = {No, No};
  T *Fij( new T(2, oo, oneBodySyms, *Cc4s::world, "Fij") );

  (*Fab)["aa"] = (*epsa)["a"];
  (*Fij)["ii"] = (*epsi)["i"];

  for (unsigned int i(0) ; i < 30 ; i++) {
    LOG(1, "Uccsd_AMPS") << "Iteration: " << i << std::endl;

    (*Rai)["bi"]  = ( - 1.0  ) * (*Fij)["ki"] * (*Tai)["bk"];
    (*Rai)["bi"] += ( + 1.0  ) * (*Fab)["bc"] * (*Tai)["ci"];
    (*Rai)["bi"] += ( - 1.0  ) * (*Tai)["cl"] * (*Viajb)["lbic"];
    (*Rai)["bi"] += ( + 0.5  ) * (*Tabij)["cblm"] * (*Vijka)["lmic"];
    (*Rai)["bi"] += ( + 0.5  ) * (*Tabij)["cdmi"] * (*Viabc)["mbcd"];
    (*Rai)["bi"] += ( - 1.0  ) * (*Tai)["bk"] * (*Tai)["dm"] * (*Vijka)["kmid"];
    (*Rai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["dm"] * (*Viabc)["mbcd"];
    (*Rai)["bi"] += ( - 0.5  ) * (*Tabij)["cblm"] * (*Tai)["fi"] * (*Vijab)["lmcf"];
    (*Rai)["bi"] += ( - 0.5  ) * (*Tabij)["cdmi"] * (*Tai)["bn"] * (*Vijab)["mncd"];
    (*Rai)["bi"] += ( + 1.0  ) * (*Tabij)["cbli"] * (*Tai)["en"] * (*Vijab)["lnce"];
    (*Rai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["bl"] * (*Tai)["en"] * (*Vijab)["lnce"];

    (*Rabij)["cdij"]  = ( + 1.0  ) * (*Vabij)["cdij"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["cm"] * (*Viajk)["mdij"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["dm"] * (*Viajk)["mcij"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Vabic)["cdie"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Vabic)["cdje"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Fij)["mi"] * (*Tabij)["cdmj"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Fij)["mj"] * (*Tabij)["cdmi"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Fab)["de"] * (*Tabij)["ecij"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Fab)["ce"] * (*Tabij)["edij"];
    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Vijkl)["mnij"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecnj"] * (*Viajb)["ndie"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ednj"] * (*Viajb)["ncie"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ecni"] * (*Viajb)["ndje"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["edni"] * (*Viajb)["ncje"];
    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["efij"] * (*Vabcd)["cdef"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijkl)["mnij"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Viajb)["ndie"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Viajb)["ncie"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Viajb)["ndje"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Viajb)["ncje"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Vabcd)["cdef"];
    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gj"] * (*Vijka)["mnig"];
    (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gi"] * (*Vijka)["mnjg"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["do"] * (*Vijka)["noie"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ednj"] * (*Tai)["co"] * (*Vijka)["noie"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Tai)["do"] * (*Vijka)["noje"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Tai)["co"] * (*Vijka)["noje"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["cdmj"] * (*Tai)["fo"] * (*Vijka)["moif"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["cdmi"] * (*Tai)["fo"] * (*Vijka)["mojf"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["gi"] * (*Viabc)["ndeg"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ednj"] * (*Tai)["gi"] * (*Viabc)["nceg"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Tai)["gj"] * (*Viabc)["ndeg"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Tai)["gj"] * (*Viabc)["nceg"];
    (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabij)["efij"] * (*Tai)["co"] * (*Viabc)["odef"];
    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["efij"] * (*Tai)["do"] * (*Viabc)["ocef"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecij"] * (*Tai)["fo"] * (*Viabc)["odef"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edij"] * (*Tai)["fo"] * (*Viabc)["ocef"];
    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["edij"] * (*Tabij)["fcop"] * (*Vijab)["opef"];
    (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabij)["ecij"] * (*Tabij)["fdop"] * (*Vijab)["opef"];
    (*Rabij)["cdij"] += ( + 0.25  ) * (*Tabij)["efij"] * (*Tabij)["cdop"] * (*Vijab)["opef"];
    (*Rabij)["cdij"] += ( - 0.5  ) * (*Tabij)["cdmi"] * (*Tabij)["fgpj"] * (*Vijab)["mpfg"];
    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tabij)["cdmj"] * (*Tabij)["fgpi"] * (*Vijab)["mpfg"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Tabij)["gcpj"] * (*Vijab)["npeg"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Tabij)["gdpj"] * (*Vijab)["npeg"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijka)["noie"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijka)["noje"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Viabc)["odef"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["do"] * (*Viabc)["ocef"];
    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tabij)["cdop"] * (*Vijab)["opef"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Tabij)["gcpj"] * (*Vijab)["npeg"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Tabij)["gdpj"] * (*Vijab)["npeg"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Tabij)["gcpi"] * (*Vijab)["npeg"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Tabij)["gdpi"] * (*Vijab)["npeg"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fo"] * (*Tabij)["cdpj"] * (*Vijab)["opef"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["fo"] * (*Tabij)["cdpi"] * (*Vijab)["opef"];
    (*Rabij)["cdij"] += ( + 0.5  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Tabij)["ghij"] * (*Vijab)["mngh"];
    (*Rabij)["cdij"] += ( - 1.0  ) * (*Tai)["dm"] * (*Tai)["fo"] * (*Tabij)["hcij"] * (*Vijab)["mofh"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["fo"] * (*Tabij)["hdij"] * (*Vijab)["mofh"];
    (*Rabij)["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Tai)["dp"] * (*Vijab)["opef"];

  }

  LOG(1,"Left EOM END") << "END" << std::endl;

}

