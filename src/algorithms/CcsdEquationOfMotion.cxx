#include <algorithms/CcsdEquationOfMotion.hpp>
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

ALGORITHM_REGISTRAR_DEFINITION(CcsdEquationOfMotion);

// TODO: Study the requirements to treat the problem with real numbers or
// complex

CcsdEquationOfMotion::CcsdEquationOfMotion(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

CcsdEquationOfMotion::~CcsdEquationOfMotion() {
}

void CcsdEquationOfMotion::run() {
  typedef CTF::Tensor<> T;

  // Get orbital energies
  T *epsi(getTensorArgument<double, T>("HoleEigenEnergies"));
  T *epsa(getTensorArgument<double, T>("ParticleEigenEnergies"));

  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int totalDimension(1 + Nv * No + No * No * Nv * Nv);
  LOG(1, "CCSD_EOM") << "Nv " << Nv << std::endl;
  LOG(1, "CCSD_EOM") << "No " << No << std::endl;
  LOG(1, "CCSD_EOM") << "Problem dimension " << totalDimension << std::endl;

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
  int oovv[] = { No, No, Nv, Nv };
  int ovoo[] = { No, Nv, No, No };
  int ovov[] = { No, Nv, No, Nv };
  int ovvv[] = { No, Nv, Nv, Nv };
  int vvov[] = { Nv, Nv, No, Nv };

  LOG(1, "CCSD_EOM") << "Antisymmetrizing Vpqrs " << std::endl;

  T *Vijab(
    new T(4, oovv, syms, *Cc4s::world, "Vijab")
  );
  (*Vijab)["ijab"] =  (*Vabij)["abij"] - (*Vabij)["abji"];

  T *Viajk(
    new T(4, ovoo, syms, *Cc4s::world, "Viajk")
  );
  (*Viajk)["iajk"] =  (*Vijka)["ijka"]  - (*Vijka)["ikja"];

  T *Viajb(
    new T(4, ovov, syms, *Cc4s::world, "Viajb")
  );
  (*Viajb)["iajb"] =  (*Vaibj)["aibj"] - (*Vaibj)["aijb"];

  T *Viabc(
    new T(4, ovvv, syms, *Cc4s::world, "Viabc")
  );
  (*Viabc)["iabc"] =  (*Vabci)["abci"] - (*Vabci)["acbi"];

  T *Vabic(
    new T(4, vvov, syms, *Cc4s::world, "Vabic")
  );
  (*Vabic)["abic"] =  (*Vabci)["abci"] - (*Vabci)["abic"];

  // Antisymmetrize integrals that are read in
  (*Vijkl)["ijkl"] -= (*Vijkl)["ijlk"];
  (*Vabcd)["abcd"] -= (*Vabcd)["abdc"];
  (*Vabij)["abij"] -= (*Vabij)["abji"];
  (*Vijka)["ijka"] -= (*Vijka)["ijak"];
  (*Vaibj)["aibj"] -= (*Vaibj)["aijb"];
  (*Vabci)["abci"] -= (*Vabci)["abic"];

  T *Tai(getTensorArgument<double, T>("CcsdSinglesAmplitudes"));
  T *Tabij(getTensorArgument<double, T>("CcsdDoublesAmplitudes"));

  CTF::Scalar<> energy(0.0);
  double energy_val(0.0);

  // Create L and R and intermediates
  int oneBodySyms[] = {NS, NS};
  int oneBodyLensL[] = {No, Nv};
  T *Lia( new T(2, oneBodyLensL, oneBodySyms, *Cc4s::world, "Lia") );
  T Lijab(false, Vijab);
  T *LHia( new T(2, oneBodyLensL, oneBodySyms, *Cc4s::world, "LHia") );
  T LHijab(false, Vijab);
  int oneBodyLensR[] = {Nv, No};
  T *Rai( new T(2, oneBodyLensR, oneBodySyms, *Cc4s::world, "Rai") );
  T Rabij(false, Vabij);
  T *HRai( new T(2, oneBodyLensR, oneBodySyms, *Cc4s::world, "HRai") );
  T HRabij(false, Vabij);

  // kinetic terms
  int kineticLensVirtual[] = {Nv, Nv};
  int kineticSyms[] = {NS, NS};
  T *Fab( new T(2, kineticLensVirtual, kineticSyms, *Cc4s::world, "Fab") );
  int kineticLensOccupied[] = {No, No};
  T *Fij( new T(2, kineticLensOccupied, kineticSyms, *Cc4s::world, "Fij") );

  // Fij and Fab are diagonal in the fock basis
  (*Fab)["aa"] = (*epsa)["a"];
  (*Fij)["ii"] = (*epsi)["i"];

  // The totalDimension should be totalDimension, but the zero-particle part is
  // zero, so we restrict the hamiltonian to the space spanned by the singles
  // and doubles excitations
  int hLens[] = {totalDimension-1, totalDimension-1};
  int hSyms[] = {NS, NS};
  T *Hpq( new CTF::Tensor<>(2, hLens, hSyms, *Cc4s::world, "Hpq") );

  int64_t *hIndices;
  double *hValues;

  for (int64_t i = 0 ; i < totalDimension-1 ; i++) {
    getCanonicalPerturbationBasis(*Lia, Lijab, i);

    /* Construct LH (One Body part) {{{ */
    (*LHia)["ja"]  = 0.0;
    (*LHia)["ja"] += ( - 1.0  ) * (*Fij)["jk"] * (*Lia)["ka"];
    (*LHia)["ja"] += ( + 1.0  ) * (*Fab)["ca"] * (*Lia)["jc"];
    (*LHia)["ja"] += ( - 1.0  ) * (*Viajb)["jcla"] * (*Lia)["lc"];
    (*LHia)["ja"] += ( - 0.5  ) * (*Viajk)["jclm"] * Lijab["mlca"];
    (*LHia)["ja"] += ( - 0.5  ) * (*Vabic)["cdma"] * Lijab["mjdc"];
    (*LHia)["ja"] += ( - 1.0  ) * (*Tai)["cl"] * (*Vijka)["ljma"] * (*Lia)["mc"];
    (*LHia)["ja"] += ( + 1.0  ) * (*Tai)["cl"] * (*Vijka)["ljmc"] * (*Lia)["ma"];
    (*LHia)["ja"] += ( - 1.0  ) * (*Tai)["cl"] * (*Viabc)["jeca"] * (*Lia)["le"];
    (*LHia)["ja"] += ( + 1.0  ) * (*Tai)["cl"] * (*Viabc)["leca"] * (*Lia)["je"];
    (*LHia)["ja"] += ( - 1.0  ) * (*Fij)["kl"] * (*Tai)["ek"] * Lijab["ljea"];
    (*LHia)["ja"] += ( + 1.0  ) * (*Fab)["cd"] * (*Tai)["dm"] * Lijab["mjca"];
    (*LHia)["ja"] += ( - 0.5  ) * (*Tai)["cl"] * (*Vijkl)["ljmn"] * Lijab["nmca"];
    (*LHia)["ja"] += ( + 1.0  ) * (*Tai)["cl"] * (*Viajb)["jenc"] * Lijab["nlea"];
    (*LHia)["ja"] += ( + 1.0  ) * (*Tai)["cl"] * (*Viajb)["lena"] * Lijab["njec"];
    (*LHia)["ja"] += ( - 1.0  ) * (*Tai)["cl"] * (*Viajb)["lenc"] * Lijab["njea"];
    (*LHia)["ja"] += ( - 0.5  ) * (*Tai)["cl"] * (*Vabcd)["efca"] * Lijab["ljfe"];
    (*LHia)["ja"] += ( + 1.0  ) * (*Tabij)["cdmn"] * (*Vijab)["njda"] * (*Lia)["mc"];
    (*LHia)["ja"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Vijab)["njcd"] * (*Lia)["ma"];
    (*LHia)["ja"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Vijab)["mnda"] * (*Lia)["jc"];
    (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Vijka)["njoa"] * Lijab["omdc"];
    (*LHia)["ja"] += ( - 1.0  ) * (*Tabij)["cdmn"] * (*Vijka)["njod"] * Lijab["omca"];
    (*LHia)["ja"] += ( - 0.25  ) * (*Tabij)["cdmn"] * (*Vijka)["mnoa"] * Lijab["ojdc"];
    (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Vijka)["mnod"] * Lijab["ojca"];
    (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Viabc)["jgda"] * Lijab["nmgc"];
    (*LHia)["ja"] += ( - 0.25  ) * (*Tabij)["cdmn"] * (*Viabc)["jgcd"] * Lijab["nmga"];
    (*LHia)["ja"] += ( - 1.0  ) * (*Tabij)["cdmn"] * (*Viabc)["ngda"] * Lijab["mjgc"];
    (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Viabc)["ngcd"] * Lijab["mjga"];
    (*LHia)["ja"] += ( - 1.0  ) * (*Tai)["cl"] * (*Tai)["en"] * (*Vijab)["njca"] * (*Lia)["le"];
    (*LHia)["ja"] += ( - 1.0  ) * (*Tai)["cl"] * (*Tai)["en"] * (*Vijab)["njec"] * (*Lia)["la"];
    (*LHia)["ja"] += ( - 1.0  ) * (*Tai)["cl"] * (*Tai)["en"] * (*Vijab)["nlea"] * (*Lia)["jc"];
    (*LHia)["ja"] += ( + 1.0  ) * (*Tai)["cl"] * (*Tai)["en"] * (*Vijka)["njoc"] * Lijab["olea"];
    (*LHia)["ja"] += ( + 0.5  ) * (*Tai)["cl"] * (*Tai)["en"] * (*Vijka)["nloa"] * Lijab["ojec"];
    (*LHia)["ja"] += ( + 1.0  ) * (*Tai)["cl"] * (*Tai)["en"] * (*Vijka)["nloe"] * Lijab["ojca"];
    (*LHia)["ja"] += ( + 0.5  ) * (*Tai)["cl"] * (*Tai)["en"] * (*Viabc)["jgec"] * Lijab["nlga"];
    (*LHia)["ja"] += ( + 1.0  ) * (*Tai)["cl"] * (*Tai)["en"] * (*Viabc)["ngca"] * Lijab["ljge"];
    (*LHia)["ja"] += ( + 1.0  ) * (*Tai)["cl"] * (*Tai)["en"] * (*Viabc)["ngec"] * Lijab["ljga"];
    (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gp"] * (*Vijab)["njga"] * Lijab["pmdc"];
    (*LHia)["ja"] += ( - 1.0  ) * (*Tabij)["cdmn"] * (*Tai)["gp"] * (*Vijab)["njgd"] * Lijab["pmca"];
    (*LHia)["ja"] += ( - 0.25  ) * (*Tabij)["cdmn"] * (*Tai)["gp"] * (*Vijab)["mnga"] * Lijab["pjdc"];
    (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gp"] * (*Vijab)["mngd"] * Lijab["pjca"];
    (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gp"] * (*Vijab)["pjda"] * Lijab["nmgc"];
    (*LHia)["ja"] += ( - 0.25  ) * (*Tabij)["cdmn"] * (*Tai)["gp"] * (*Vijab)["pjcd"] * Lijab["nmga"];
    (*LHia)["ja"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gp"] * (*Vijab)["pjgd"] * Lijab["nmca"];
    (*LHia)["ja"] += ( - 1.0  ) * (*Tabij)["cdmn"] * (*Tai)["gp"] * (*Vijab)["pnda"] * Lijab["mjgc"];
    (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gp"] * (*Vijab)["pncd"] * Lijab["mjga"];
    (*LHia)["ja"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gp"] * (*Vijab)["pnga"] * Lijab["mjdc"];
    (*LHia)["ja"] += ( + 1.0  ) * (*Tabij)["cdmn"] * (*Tai)["gp"] * (*Vijab)["pngd"] * Lijab["mjca"];
    (*LHia)["ja"] += ( + 0.5  ) * (*Tai)["cl"] * (*Tai)["en"] * (*Tai)["gp"] * (*Vijab)["pjec"] * Lijab["nlga"];
    (*LHia)["ja"] += ( + 0.5  ) * (*Tai)["cl"] * (*Tai)["en"] * (*Tai)["gp"] * (*Vijab)["pnca"] * Lijab["ljge"];
    (*LHia)["ja"] += ( - 1.0  ) * (*Tai)["cl"] * (*Tai)["en"] * (*Tai)["gp"] * (*Vijab)["pngc"] * Lijab["ljea"];
    /* Construct LH (One Body part) }}} */

    /* Construct LH (Two Body part) {{{ */
    LHijab["klab"]  = 0.0;
    LHijab["klab"] += ( - 1.0  ) * (*Vijka)["klmb"] * (*Lia)["ma"];
    LHijab["klab"] += ( + 1.0  ) * (*Vijka)["klma"] * (*Lia)["mb"];
    LHijab["klab"] += ( + 1.0  ) * (*Viabc)["keab"] * (*Lia)["le"];
    LHijab["klab"] += ( - 1.0  ) * (*Viabc)["leab"] * (*Lia)["ke"];
    LHijab["klab"] += ( - 1.0  ) * (*Fij)["km"] * Lijab["mlab"];
    LHijab["klab"] += ( + 1.0  ) * (*Fij)["lm"] * Lijab["mkab"];
    LHijab["klab"] += ( - 1.0  ) * (*Fab)["eb"] * Lijab["klea"];
    LHijab["klab"] += ( + 1.0  ) * (*Fab)["ea"] * Lijab["kleb"];
    LHijab["klab"] += ( - 0.5  ) * (*Vijkl)["klmn"] * Lijab["nmab"];
    LHijab["klab"] += ( + 1.0  ) * (*Viajb)["kenb"] * Lijab["nlea"];
    LHijab["klab"] += ( - 1.0  ) * (*Viajb)["kena"] * Lijab["nleb"];
    LHijab["klab"] += ( - 1.0  ) * (*Viajb)["lenb"] * Lijab["nkea"];
    LHijab["klab"] += ( + 1.0  ) * (*Viajb)["lena"] * Lijab["nkeb"];
    LHijab["klab"] += ( - 0.5  ) * (*Vabcd)["efab"] * Lijab["klfe"];
    LHijab["klab"] += ( - 1.0  ) * (*Tai)["en"] * (*Vijab)["kleb"] * (*Lia)["na"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Vijab)["klea"] * (*Lia)["nb"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Vijab)["nkab"] * (*Lia)["le"];
    LHijab["klab"] += ( - 1.0  ) * (*Tai)["en"] * (*Vijab)["nlab"] * (*Lia)["ke"];
    LHijab["klab"] += ( - 1.0  ) * (*Tai)["en"] * (*Vijab)["nkeb"] * (*Lia)["la"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Vijab)["nkea"] * (*Lia)["lb"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Vijab)["nleb"] * (*Lia)["ka"];
    LHijab["klab"] += ( - 1.0  ) * (*Tai)["en"] * (*Vijab)["nlea"] * (*Lia)["kb"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Vijka)["kloe"] * Lijab["onab"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Vijka)["nkob"] * Lijab["olea"];
    LHijab["klab"] += ( - 1.0  ) * (*Tai)["en"] * (*Vijka)["nkoa"] * Lijab["oleb"];
    LHijab["klab"] += ( - 1.0  ) * (*Tai)["en"] * (*Vijka)["nlob"] * Lijab["okea"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Vijka)["nloa"] * Lijab["okeb"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Vijka)["nkoe"] * Lijab["olab"];
    LHijab["klab"] += ( - 1.0  ) * (*Tai)["en"] * (*Vijka)["nloe"] * Lijab["okab"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Viabc)["kgeb"] * Lijab["nlga"];
    LHijab["klab"] += ( - 1.0  ) * (*Tai)["en"] * (*Viabc)["kgea"] * Lijab["nlgb"];
    LHijab["klab"] += ( - 1.0  ) * (*Tai)["en"] * (*Viabc)["lgeb"] * Lijab["nkga"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Viabc)["lgea"] * Lijab["nkgb"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Viabc)["ngab"] * Lijab["klge"];
    LHijab["klab"] += ( - 1.0  ) * (*Tai)["en"] * (*Viabc)["ngeb"] * Lijab["klga"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Viabc)["ngea"] * Lijab["klgb"];
    LHijab["klab"] += ( + 0.5  ) * (*Tabij)["efop"] * (*Vijab)["klfb"] * Lijab["poea"];
    LHijab["klab"] += ( - 0.5  ) * (*Tabij)["efop"] * (*Vijab)["klfa"] * Lijab["poeb"];
    LHijab["klab"] += ( - 0.25  ) * (*Tabij)["efop"] * (*Vijab)["klef"] * Lijab["poab"];
    LHijab["klab"] += ( - 0.5  ) * (*Tabij)["efop"] * (*Vijab)["pkab"] * Lijab["olfe"];
    LHijab["klab"] += ( + 0.5  ) * (*Tabij)["efop"] * (*Vijab)["plab"] * Lijab["okfe"];
    LHijab["klab"] += ( - 1.0  ) * (*Tabij)["efop"] * (*Vijab)["pkfb"] * Lijab["olea"];
    LHijab["klab"] += ( + 1.0  ) * (*Tabij)["efop"] * (*Vijab)["pkfa"] * Lijab["oleb"];
    LHijab["klab"] += ( + 1.0  ) * (*Tabij)["efop"] * (*Vijab)["plfb"] * Lijab["okea"];
    LHijab["klab"] += ( - 1.0  ) * (*Tabij)["efop"] * (*Vijab)["plfa"] * Lijab["okeb"];
    LHijab["klab"] += ( + 0.5  ) * (*Tabij)["efop"] * (*Vijab)["pkef"] * Lijab["olab"];
    LHijab["klab"] += ( - 0.5  ) * (*Tabij)["efop"] * (*Vijab)["plef"] * Lijab["okab"];
    LHijab["klab"] += ( - 0.25  ) * (*Tabij)["efop"] * (*Vijab)["opab"] * Lijab["klfe"];
    LHijab["klab"] += ( - 0.5  ) * (*Tabij)["efop"] * (*Vijab)["opfb"] * Lijab["klea"];
    LHijab["klab"] += ( + 0.5  ) * (*Tabij)["efop"] * (*Vijab)["opfa"] * Lijab["kleb"];
    LHijab["klab"] += ( + 0.5  ) * (*Tai)["en"] * (*Tai)["gp"] * (*Vijab)["klge"] * Lijab["pnab"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Tai)["gp"] * (*Vijab)["pkeb"] * Lijab["nlga"];
    LHijab["klab"] += ( - 1.0  ) * (*Tai)["en"] * (*Tai)["gp"] * (*Vijab)["pkea"] * Lijab["nlgb"];
    LHijab["klab"] += ( - 1.0  ) * (*Tai)["en"] * (*Tai)["gp"] * (*Vijab)["pleb"] * Lijab["nkga"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Tai)["gp"] * (*Vijab)["plea"] * Lijab["nkgb"];
    LHijab["klab"] += ( - 1.0  ) * (*Tai)["en"] * (*Tai)["gp"] * (*Vijab)["pkge"] * Lijab["nlab"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Tai)["gp"] * (*Vijab)["plge"] * Lijab["nkab"];
    LHijab["klab"] += ( + 0.5  ) * (*Tai)["en"] * (*Tai)["gp"] * (*Vijab)["pnab"] * Lijab["klge"];
    LHijab["klab"] += ( + 1.0  ) * (*Tai)["en"] * (*Tai)["gp"] * (*Vijab)["pngb"] * Lijab["klea"];
    LHijab["klab"] += ( - 1.0  ) * (*Tai)["en"] * (*Tai)["gp"] * (*Vijab)["pnga"] * Lijab["kleb"];
    /* }}} Construct LH (Two Body part) */

    for (int64_t j = 0 ; j < totalDimension-1; j++) {
      getCanonicalPerturbationBasis(*Rai, Rabij, j);

      /* Construct HR (One Body part) {{{ */
      (*HRai)["bi"]  = 0.0;
      (*HRai)["bi"] += ( - 1.0  ) * (*Fij)["ki"] * (*Rai)["bk"];
      (*HRai)["bi"] += ( + 1.0  ) * (*Fab)["bc"] * (*Rai)["ci"];
      (*HRai)["bi"] += ( - 1.0  ) * (*Viajb)["kbid"] * (*Rai)["dk"];
      (*HRai)["bi"] += ( + 0.5  ) * (*Vijka)["klie"] * Rabij["ebkl"];
      (*HRai)["bi"] += ( + 0.5  ) * (*Viabc)["kbde"] * Rabij["deki"];
      (*HRai)["bi"] += ( - 1.0  ) * (*Tai)["bk"] * (*Vijka)["klie"] * (*Rai)["el"];
      (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["cl"] * (*Vijka)["lmic"] * (*Rai)["bm"];
      (*HRai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Viabc)["lbce"] * (*Rai)["el"];
      (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["cl"] * (*Viabc)["lbce"] * (*Rai)["ei"];
      (*HRai)["bi"] += ( + 0.5  ) * (*Tai)["ci"] * (*Vijab)["lmcf"] * Rabij["fblm"];
      (*HRai)["bi"] += ( + 0.5  ) * (*Tai)["bk"] * (*Vijab)["klef"] * Rabij["efli"];
      (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["cl"] * (*Vijab)["lmcf"] * Rabij["fbmi"];
      (*HRai)["bi"] += ( + 1.0  ) * (*Tabij)["cbli"] * (*Vijab)["lmcf"] * (*Rai)["fm"];
      (*HRai)["bi"] += ( - 0.5  ) * (*Tabij)["cdmi"] * (*Vijab)["mncd"] * (*Rai)["bn"];
      (*HRai)["bi"] += ( - 0.5  ) * (*Tabij)["cblm"] * (*Vijab)["lmcf"] * (*Rai)["fi"];
      (*HRai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["bl"] * (*Vijab)["lmcf"] * (*Rai)["fm"];
      (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["ci"] * (*Tai)["dm"] * (*Vijab)["mncd"] * (*Rai)["bn"];
      (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["bk"] * (*Tai)["dm"] * (*Vijab)["kmdf"] * (*Rai)["fi"];
      /* }}} Construct HR (One Body part) */

      /* Construct HR (One Body part) {{{ */
      HRabij["cdij"]  = 0.0;
      HRabij["cdij"] += ( - 1.0  ) * (*Viajk)["mdij"] * (*Rai)["cm"];
      HRabij["cdij"] += ( + 1.0  ) * (*Viajk)["mcij"] * (*Rai)["dm"];
      HRabij["cdij"] += ( + 1.0  ) * (*Vabic)["cdie"] * (*Rai)["ej"];
      HRabij["cdij"] += ( - 1.0  ) * (*Vabic)["cdje"] * (*Rai)["ei"];
      HRabij["cdij"] += ( - 1.0  ) * (*Fij)["mi"] * Rabij["cdmj"];
      HRabij["cdij"] += ( + 1.0  ) * (*Fij)["mj"] * Rabij["cdmi"];
      HRabij["cdij"] += ( - 1.0  ) * (*Fab)["de"] * Rabij["ecij"];
      HRabij["cdij"] += ( + 1.0  ) * (*Fab)["ce"] * Rabij["edij"];
      HRabij["cdij"] += ( + 0.5  ) * (*Vijkl)["mnij"] * Rabij["cdmn"];
      HRabij["cdij"] += ( + 1.0  ) * (*Viajb)["mdif"] * Rabij["fcmj"];
      HRabij["cdij"] += ( - 1.0  ) * (*Viajb)["mcif"] * Rabij["fdmj"];
      HRabij["cdij"] += ( - 1.0  ) * (*Viajb)["mdjf"] * Rabij["fcmi"];
      HRabij["cdij"] += ( + 1.0  ) * (*Viajb)["mcjf"] * Rabij["fdmi"];
      HRabij["cdij"] += ( + 0.5  ) * (*Vabcd)["cdef"] * Rabij["efij"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["dm"] * (*Vijkl)["mnij"] * (*Rai)["cn"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["cm"] * (*Vijkl)["mnij"] * (*Rai)["dn"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Viajb)["ndie"] * (*Rai)["cn"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Viajb)["ncie"] * (*Rai)["dn"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Viajb)["ndje"] * (*Rai)["cn"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Viajb)["ncje"] * (*Rai)["dn"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["cm"] * (*Viajb)["mdif"] * (*Rai)["fj"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["dm"] * (*Viajb)["mcif"] * (*Rai)["fj"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["cm"] * (*Viajb)["mdjf"] * (*Rai)["fi"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["dm"] * (*Viajb)["mcjf"] * (*Rai)["fi"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Vabcd)["cdef"] * (*Rai)["fj"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Vabcd)["cdef"] * (*Rai)["fi"];
      HRabij["cdij"] += ( + 0.5  ) * (*Tai)["ej"] * (*Vijka)["noie"] * Rabij["cdno"];
      HRabij["cdij"] += ( - 0.5  ) * (*Tai)["ei"] * (*Vijka)["noje"] * Rabij["cdno"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["dm"] * (*Vijka)["mnig"] * Rabij["gcnj"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["cm"] * (*Vijka)["mnig"] * Rabij["gdnj"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["dm"] * (*Vijka)["mnjg"] * Rabij["gcni"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["cm"] * (*Vijka)["mnjg"] * Rabij["gdni"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["en"] * (*Vijka)["noie"] * Rabij["cdoj"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["en"] * (*Vijka)["noje"] * Rabij["cdoi"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Viabc)["ndeg"] * Rabij["gcnj"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Viabc)["nceg"] * Rabij["gdnj"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Viabc)["ndeg"] * Rabij["gcni"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Viabc)["nceg"] * Rabij["gdni"];
      HRabij["cdij"] += ( - 0.5  ) * (*Tai)["cm"] * (*Viabc)["mdfg"] * Rabij["fgij"];
      HRabij["cdij"] += ( + 0.5  ) * (*Tai)["dm"] * (*Viabc)["mcfg"] * Rabij["fgij"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["en"] * (*Viabc)["ndeg"] * Rabij["gcij"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["en"] * (*Viabc)["nceg"] * Rabij["gdij"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["cdmj"] * (*Vijka)["mnig"] * (*Rai)["gn"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["cdmi"] * (*Vijka)["mnjg"] * (*Rai)["gn"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["ednj"] * (*Vijka)["noie"] * (*Rai)["co"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["ecnj"] * (*Vijka)["noie"] * (*Rai)["do"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Vijka)["noje"] * (*Rai)["co"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Vijka)["noje"] * (*Rai)["do"];
      HRabij["cdij"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Vijka)["mnig"] * (*Rai)["gj"];
      HRabij["cdij"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Vijka)["mnjg"] * (*Rai)["gi"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["ecij"] * (*Viabc)["ndeg"] * (*Rai)["gn"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["edij"] * (*Viabc)["nceg"] * (*Rai)["gn"];
      HRabij["cdij"] += ( - 0.5  ) * (*Tabij)["efij"] * (*Viabc)["odef"] * (*Rai)["co"];
      HRabij["cdij"] += ( + 0.5  ) * (*Tabij)["efij"] * (*Viabc)["ocef"] * (*Rai)["do"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Viabc)["ndeg"] * (*Rai)["gj"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Viabc)["nceg"] * (*Rai)["gj"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["ecnj"] * (*Viabc)["ndeg"] * (*Rai)["gi"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["ednj"] * (*Viabc)["nceg"] * (*Rai)["gi"];
      HRabij["cdij"] += ( + 0.5  ) * (*Tabij)["edij"] * (*Vijab)["noeh"] * Rabij["hcno"];
      HRabij["cdij"] += ( - 0.5  ) * (*Tabij)["ecij"] * (*Vijab)["noeh"] * Rabij["hdno"];
      HRabij["cdij"] += ( + 0.25  ) * (*Tabij)["efij"] * (*Vijab)["opef"] * Rabij["cdop"];
      HRabij["cdij"] += ( - 0.5  ) * (*Tabij)["cdmi"] * (*Vijab)["mngh"] * Rabij["ghnj"];
      HRabij["cdij"] += ( + 0.5  ) * (*Tabij)["cdmj"] * (*Vijab)["mngh"] * Rabij["ghni"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["edni"] * (*Vijab)["noeh"] * Rabij["hcoj"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["ecni"] * (*Vijab)["noeh"] * Rabij["hdoj"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["ednj"] * (*Vijab)["noeh"] * Rabij["hcoi"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["ecnj"] * (*Vijab)["noeh"] * Rabij["hdoi"];
      HRabij["cdij"] += ( - 0.5  ) * (*Tabij)["efoi"] * (*Vijab)["opef"] * Rabij["cdpj"];
      HRabij["cdij"] += ( + 0.5  ) * (*Tabij)["efoj"] * (*Vijab)["opef"] * Rabij["cdpi"];
      HRabij["cdij"] += ( + 0.25  ) * (*Tabij)["cdmn"] * (*Vijab)["mngh"] * Rabij["ghij"];
      HRabij["cdij"] += ( + 0.5  ) * (*Tabij)["edno"] * (*Vijab)["noeh"] * Rabij["hcij"];
      HRabij["cdij"] += ( - 0.5  ) * (*Tabij)["ecno"] * (*Vijab)["noeh"] * Rabij["hdij"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Vijka)["noie"] * (*Rai)["co"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Vijka)["noie"] * (*Rai)["do"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Vijka)["noje"] * (*Rai)["co"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Vijka)["noje"] * (*Rai)["do"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijka)["mnig"] * (*Rai)["gj"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijka)["mnjg"] * (*Rai)["gi"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Viabc)["odef"] * (*Rai)["co"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Viabc)["ocef"] * (*Rai)["do"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Viabc)["ndeg"] * (*Rai)["gj"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Viabc)["nceg"] * (*Rai)["gj"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Viabc)["ndeg"] * (*Rai)["gi"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Viabc)["nceg"] * (*Rai)["gi"];
      HRabij["cdij"] += ( + 0.5  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Vijab)["opef"] * Rabij["cdop"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Vijab)["noeh"] * Rabij["hcoj"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Vijab)["noeh"] * Rabij["hdoj"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Vijab)["noeh"] * Rabij["hcoi"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Vijab)["noeh"] * Rabij["hdoi"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fo"] * (*Vijab)["opef"] * Rabij["cdpj"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["fo"] * (*Vijab)["opef"] * Rabij["cdpi"];
      HRabij["cdij"] += ( + 0.5  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijab)["mngh"] * Rabij["ghij"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["dm"] * (*Tai)["fo"] * (*Vijab)["mofh"] * Rabij["hcij"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["fo"] * (*Vijab)["mofh"] * Rabij["hdij"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["cdmj"] * (*Tai)["fi"] * (*Vijab)["mofh"] * (*Rai)["ho"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["cdmi"] * (*Tai)["fj"] * (*Vijab)["mofh"] * (*Rai)["ho"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["ednj"] * (*Tai)["gi"] * (*Vijab)["npeg"] * (*Rai)["cp"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["gi"] * (*Vijab)["npeg"] * (*Rai)["dp"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["edni"] * (*Tai)["gj"] * (*Vijab)["npeg"] * (*Rai)["cp"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["ecni"] * (*Tai)["gj"] * (*Vijab)["npeg"] * (*Rai)["dp"];
      HRabij["cdij"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gi"] * (*Vijab)["mngh"] * (*Rai)["hj"];
      HRabij["cdij"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gj"] * (*Vijab)["mngh"] * (*Rai)["hi"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["ecij"] * (*Tai)["dn"] * (*Vijab)["noeh"] * (*Rai)["ho"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["edij"] * (*Tai)["cn"] * (*Vijab)["noeh"] * (*Rai)["ho"];
      HRabij["cdij"] += ( - 0.5  ) * (*Tabij)["efij"] * (*Tai)["do"] * (*Vijab)["opef"] * (*Rai)["cp"];
      HRabij["cdij"] += ( + 0.5  ) * (*Tabij)["efij"] * (*Tai)["co"] * (*Vijab)["opef"] * (*Rai)["dp"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["edij"] * (*Tai)["fo"] * (*Vijab)["opef"] * (*Rai)["cp"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["ecij"] * (*Tai)["fo"] * (*Vijab)["opef"] * (*Rai)["dp"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["ecni"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hj"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["edni"] * (*Tai)["co"] * (*Vijab)["noeh"] * (*Rai)["hj"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hi"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["ednj"] * (*Tai)["co"] * (*Vijab)["noeh"] * (*Rai)["hi"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tabij)["cdmi"] * (*Tai)["fo"] * (*Vijab)["mofh"] * (*Rai)["hj"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tabij)["cdmj"] * (*Tai)["fo"] * (*Vijab)["mofh"] * (*Rai)["hi"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["do"] * (*Vijab)["opef"] * (*Rai)["cp"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Vijab)["opef"] * (*Rai)["dp"];
      HRabij["cdij"] += ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hj"];
      HRabij["cdij"] += ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hi"];
      /* }}} Construct HR (Two Body part) */


      energy[""]  = (*LHia)["ia"] * (*HRai)["ai"];
      energy[""] += LHijab["ijab"] * HRabij["abij"];
      energy_val = energy.get_val();
      int Hpq_data_size;

      if (Hpq->wrld->rank == 0) {
        Hpq_data_size = 1;
        hValues = (double*) malloc(1);
        hIndices = (int64_t*) malloc(1);
        hValues[0] = energy_val;
        hIndices[0] = i + j * (totalDimension - 1);
      } else {
        Hpq_data_size = 0;
        hValues = (double*) malloc(0);
        hIndices = (int64_t*) malloc(0);
      }

      (*Hpq).write(Hpq_data_size, hIndices, hValues);
      LOG(1, "CCSD_EOM") << "< " << i << " |H| " <<  j << " >"  << " = "
                        << energy_val << std::endl;
    }
  }

  allocatedTensorArgument("SimlarityTransformedHamiltonianSD", Hpq);

}

template <typename F>
void CcsdEquationOfMotion::getCanonicalPerturbationBasis(
    CTF::Tensor<F> &Tai, CTF::Tensor<F> &Tabij, int64_t i) {
  int oneBodyLength(Tai.lens[0] * Tai.lens[1]);
  int twoBodyLength(
      Tabij.lens[0] * Tabij.lens[1] * Tabij.lens[2] *  Tabij.lens[3]
  );

  Tabij["abij"] = 0;
  Tai["ai"] = 0;
  int64_t *indices;
  F *values;
  int arrayCount;

  if (Tabij.wrld->rank == 0) {
    arrayCount = 1;
    values = (F*) malloc(arrayCount);
    indices = (int64_t*) malloc(arrayCount);
    indices[0] = i + 1 <= oneBodyLength ? i : i - oneBodyLength;
    values[0] = 1.0;
  } else {
    arrayCount = 0;
    values = (F*) malloc(arrayCount);
    indices = (int64_t*) malloc(arrayCount);
  }

  if (i+1 <= oneBodyLength) { // One body regime
    Tai.write(arrayCount, indices, values);
  } else { // Two body regime
    Tabij.write(arrayCount, indices, values);
  }
  //Tai.print();
  //Tabij.print();

}

// instantiate
template
void CcsdEquationOfMotion::getCanonicalPerturbationBasis(
    CTF::Tensor<double> &Tai, CTF::Tensor<double> &Tabij, int64_t i);
