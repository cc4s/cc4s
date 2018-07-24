#include <algorithms/CcsdEquationOfMotionDavidson.hpp>

#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <math/EigenSystemDavidson.hpp>
#include <math/MathFunctions.hpp>
#include <math/FockVector.hpp>
#include <math/ComplexTensor.hpp>
#include <math/RandomTensor.hpp>
#include <util/MpiCommunicator.hpp>
#include <util/Log.hpp>
#include <util/TensorIo.hpp>
#include <util/Exception.hpp>
#include <util/RangeParser.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <util/SharedPointer.hpp>

#include <algorithm>
#include <utility>
#include <limits>

using namespace cc4s;
using namespace tcc;

ALGORITHM_REGISTRAR_DEFINITION(CcsdEquationOfMotionDavidson);

CcsdEquationOfMotionDavidson::CcsdEquationOfMotionDavidson(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}
CcsdEquationOfMotionDavidson::~CcsdEquationOfMotionDavidson() {}

void CcsdEquationOfMotionDavidson::run() {

  // Get copy of couloumb integrals
  CTF::Tensor<> *Vijkl(
      getTensorArgument<double, CTF::Tensor<> >("HHHHCoulombIntegrals"));
  CTF::Tensor<> *Vabcd(
      getTensorArgument<double, CTF::Tensor<> >("PPPPCoulombIntegrals"));
  CTF::Tensor<> *Vijka(
      getTensorArgument<double, CTF::Tensor<> >("HHHPCoulombIntegrals"));
  CTF::Tensor<> *Vijab(
      getTensorArgument<double, CTF::Tensor<> >("HHPPCoulombIntegrals"));
  CTF::Tensor<> *Viajk(
      getTensorArgument<double, CTF::Tensor<> >("HPHHCoulombIntegrals"));
  CTF::Tensor<> *Viajb(
      getTensorArgument<double, CTF::Tensor<> >("HPHPCoulombIntegrals"));
  CTF::Tensor<> *Viabc(
      getTensorArgument<double, CTF::Tensor<> >("HPPPCoulombIntegrals"));
  CTF::Tensor<> *Vabic(
      getTensorArgument<double, CTF::Tensor<> >("PPHPCoulombIntegrals"));
  CTF::Tensor<> *Vabci(
      getTensorArgument<double, CTF::Tensor<> >("PPPHCoulombIntegrals"));
  CTF::Tensor<> *Vaibc(
      getTensorArgument<double, CTF::Tensor<> >("PHPPCoulombIntegrals"));
  CTF::Tensor<> *Vaibj(
      getTensorArgument<double, CTF::Tensor<> >("PHPHCoulombIntegrals"));
  CTF::Tensor<> *Viabj(
      getTensorArgument<double, CTF::Tensor<> >("HPPHCoulombIntegrals"));
  CTF::Tensor<> *Vijak(
      getTensorArgument<double, CTF::Tensor<> >("HHPHCoulombIntegrals"));
  CTF::Tensor<> *Vaijb(
      getTensorArgument<double, CTF::Tensor<> >("PHHPCoulombIntegrals"));

  //CTF::Tensor<> *Vabij(
      //getTensorArgument<double, CTF::Tensor<>>("PPHHCoulombIntegrals"));

  // Get orbital energies
  CTF::Tensor<> *epsi(
      getTensorArgument<double, CTF::Tensor<> >("HoleEigenEnergies"));
  CTF::Tensor<> *epsa(
      getTensorArgument<double, CTF::Tensor<> >("ParticleEigenEnergies"));
  int Nv(epsa->lens[0]), No(epsi->lens[0]);

  // HF terms
  int vv[] = {Nv, Nv};
  int oo[] = {No, No};
  int kineticSyms[] = {NS, NS};
  CTF::Tensor<> *Fab(
    new CTF::Tensor<>(2, vv, kineticSyms, *Cc4s::world, "Fab")
  );
  CTF::Tensor<> *Fij(
    new CTF::Tensor<>(2, oo, kineticSyms, *Cc4s::world, "Fij")
  );
  CTF::Tensor<> *Fia;

  if (
    isArgumentGiven("HPFockMatrix") &&
    isArgumentGiven("HHFockMatrix") &&
    isArgumentGiven("PPFockMatrix")
  ) {
    LOG(0, "CcsdEomDavid") << "Using non-canonical orbitals" << std::endl;
    Fia = getTensorArgument<double, CTF::Tensor<> >("HPFockMatrix");
    Fab = getTensorArgument<double, CTF::Tensor<> >("PPFockMatrix");
    Fij = getTensorArgument<double, CTF::Tensor<> >("HHFockMatrix");
  } else {
    LOG(0, "CcsdEomDavid") << "Using canonical orbitals" << std::endl;
    Fia = nullptr;
    (*Fab)["aa"] = (*epsa)["a"];
    (*Fij)["ii"] = (*epsi)["i"];
  }



  int syms2[] = {NS, NS};
  int syms4[] = {NS, NS, NS, NS};
  int vo[] = {Nv,No};
  int vvoo[] = {Nv,Nv,No,No};
  // We initialize the T amplitudes here so that it is not necessary
  // to do a ccsd calculation before to do the CISD calculation.
  CTF::Tensor<> Tai(2, vo, syms2, *Cc4s::world, "Tai");
  CTF::Tensor<> Tabij(4, vvoo, syms4, *Cc4s::world, "Tabij");
  if (getIntegerArgument("CISD", 0) == 1) {
    LOG(0, "CcsdEomDavid") << "Calculating CISD" << std::endl;
    Tai["ai"] = 0.0;
    Tabij["abij"] = 0.0;
  } else {
    // Get the Uccsd amplitudes from the input file
    Tai["ai"] =
    (*getTensorArgument<double, CTF::Tensor<> >("SinglesAmplitudes"))["ai"];
    Tabij["abij"] =
    (*getTensorArgument<double, CTF::Tensor<> >("DoublesAmplitudes"))["abij"];
  }

  if (getIntegerArgument("printTensors", 0) == 1) {
    TensorIo::writeText<>(
      "Tai.tensor", Tai, "ij", "", " "
    );
    TensorIo::writeText<>(
      "Tabij.tensor", Tabij, "ijkl", "", " "
    );
  }

  CcsdSimilarityTransformedHamiltonian<double> H(
    &Tai, &Tabij, Fij, Fab, Fia,
    Vabcd, Viajb, Vijab, Vijkl, Vijka, Viabc, Viajk, Vabic,
    Vaibc, Vaibj, Viabj, Vijak, Vaijb, Vabci
  );

  unsigned int maxIterations(getIntegerArgument("maxIterations", 32));
  unsigned int minIterations(getIntegerArgument("minIterations", 1));
  bool intermediates(getIntegerArgument("intermediates", 1));
  bool singleParticleOnly(getIntegerArgument("singleParticleOnly", 0));
  H.buildIntermediates(intermediates, singleParticleOnly);

  CcsdPreConditioner<double> P(
    Tai, Tabij, *Fij, *Fab, *Vabcd, *Viajb, *Vijab, *Vijkl, this
  );
  allocatedTensorArgument(
    "SinglesHamiltonianDiagonal",
    new CTF::Tensor<>(*P.getDiagonalH().get(0))
  );
  allocatedTensorArgument(
    "DoublesHamiltonianDiagonal",
    new CTF::Tensor<>(*P.getDiagonalH().get(1))
  );

  // Davidson solver
  int eigenStates(getIntegerArgument("eigenstates", 1));
  LOG(0, "CcsdEomDavid") << "Max iterations " << maxIterations << std::endl;
  double ediff(getRealArgument("ediff", 1e-4));
  LOG(0, "CcsdEomDavid") << "ediff " << ediff << std::endl;
  LOG(0, "CcsdEomDavid") << "Computing " << eigenStates << " eigen states"
                              << std::endl;
  std::vector<int> refreshIterations(
    RangeParser(getTextArgument("refreshIterations", "")).getRange()
  );
  EigenSystemDavidson<FockVector<double>> eigenSystem(
    H, eigenStates, P, ediff,
    No*Nv + (No*(No - 1)/2 ) * (Nv * (Nv - 1)/2),
    maxIterations, minIterations, false, refreshIterations
  );

  std::vector<complex> eigenValues(eigenSystem.getEigenValues());
  int eigenCounter(0);
  for (auto &ev: eigenValues) {
    eigenCounter++;
    LOG(0, "CcsdEomDavid") << eigenCounter << ". Eigenvalue=" << ev << std::endl;
  }
  if (getIntegerArgument("printTensors", 0) == 1) {
    eigenCounter = 0;
    for (auto &eigenState: eigenSystem.getRightEigenVectors()) {
      eigenCounter++;
      LOG(0, "CcsdEomDavid") << eigenCounter << ". Eigenstate=" << std::endl;
      eigenState.get(1)->print(stdout, -1e100);
      eigenState.get(2)->print(stdout, -1e100);
    }
  }
}


template <typename F>
CcsdSimilarityTransformedHamiltonian<F>::CcsdSimilarityTransformedHamiltonian(
  CTF::Tensor<F> *Tai_,
  CTF::Tensor<F> *Tabij_,
  CTF::Tensor<F> *Fij_,
  CTF::Tensor<F> *Fab_,
  CTF::Tensor<F> *Fia_,
  CTF::Tensor<F> *Vabcd_,
  CTF::Tensor<F> *Viajb_,
  CTF::Tensor<F> *Vijab_,
  CTF::Tensor<F> *Vijkl_,
  CTF::Tensor<F> *Vijka_,
  CTF::Tensor<F> *Viabc_,
  CTF::Tensor<F> *Viajk_,
  CTF::Tensor<F> *Vabic_,
  CTF::Tensor<F> *Vaibc_,
  CTF::Tensor<F> *Vaibj_,
  CTF::Tensor<F> *Viabj_,
  CTF::Tensor<F> *Vijak_,
  CTF::Tensor<F> *Vaijb_,
  CTF::Tensor<F> *Vabci_
):
  Tai(Tai_),
  Tabij(Tabij_),
  Fij(Fij_),
  Fab(Fab_),
  Fia(Fia_),
  Vabcd(Vabcd_),
  Viajb(Viajb_),
  Vijab(Vijab_),
  Vijkl(Vijkl_),
  Vijka(Vijka_),
  Viabc(Viabc_),
  Viajk(Viajk_),
  Vabic(Vabic_),
  Vaibc(Vaibc_),
  Vaibj(Vaibj_),
  Viabj(Viabj_),
  Vijak(Vijak_),
  Vaijb(Vaijb_),
  Vabci(Vabci_)
{
}

template <typename F>
FockVector<F> CcsdSimilarityTransformedHamiltonian<F>::leftApply(
  FockVector<F> &L
) {
  FockVector<F> LH(L);
  // get pointers to the component tensors
  PTR(CTF::Tensor<F>) Lia( L.get(0) );
  PTR(CTF::Tensor<F>) Lijab( L.get(1) );
  PTR(CTF::Tensor<F>) LHia( LH.get(0) );
  PTR(CTF::Tensor<F>) LHijab( LH.get(1) );

  // Contruct HR (one body part)
  (*LHia)["ja"]  = 0.0;
  (*LHia)["ja"] += ( - 1.0  ) * (*Fij)["jk"] * (*Lia)["ka"];
  (*LHia)["ja"] += ( + 1.0  ) * (*Fab)["ca"] * (*Lia)["jc"];
  (*LHia)["ja"] += ( - 1.0  ) * (*Viajb)["jcla"] * (*Lia)["lc"];
  (*LHia)["ja"] += ( - 0.5  ) * (*Viajk)["jclm"] * (*Lijab)["mlca"];
  (*LHia)["ja"] += ( - 0.5  ) * (*Vabic)["cdma"] * (*Lijab)["mjdc"];
  (*LHia)["ja"] += ( + 1.0  ) * (*Tabij)["cdmn"] * (*Vijab)["njda"] * (*Lia)["mc"];
  (*LHia)["ja"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Vijab)["njcd"] * (*Lia)["ma"];
  (*LHia)["ja"] += ( + 0.5  ) * (*Tabij)["cdmn"] * (*Vijab)["mnda"] * (*Lia)["jc"];
  (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Vijka)["njoa"] * (*Lijab)["omdc"];
  (*LHia)["ja"] += ( - 1.0  ) * (*Tabij)["cdmn"] * (*Vijka)["njod"] * (*Lijab)["omca"];
  (*LHia)["ja"] += ( - 0.25  ) * (*Tabij)["cdmn"] * (*Vijka)["mnoa"] * (*Lijab)["ojdc"];
  (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Vijka)["mnod"] * (*Lijab)["ojca"];
  (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Viabc)["jgda"] * (*Lijab)["nmgc"];
  (*LHia)["ja"] += ( - 0.25  ) * (*Tabij)["cdmn"] * (*Viabc)["jgcd"] * (*Lijab)["nmga"];
  (*LHia)["ja"] += ( - 1.0  ) * (*Tabij)["cdmn"] * (*Viabc)["ngda"] * (*Lijab)["mjgc"];
  (*LHia)["ja"] += ( - 0.5  ) * (*Tabij)["cdmn"] * (*Viabc)["ngcd"] * (*Lijab)["mjga"];

  // Contruct HR (two body part)
  (*LHijab)["klab"]  = 0.0;
  (*LHijab)["klab"] += ( - 1.0  ) * (*Vijka)["klmb"] * (*Lia)["ma"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Vijka)["klma"] * (*Lia)["mb"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Viabc)["keab"] * (*Lia)["le"];
  (*LHijab)["klab"] += ( - 1.0  ) * (*Viabc)["leab"] * (*Lia)["ke"];
  (*LHijab)["klab"] += ( - 1.0  ) * (*Fij)["km"] * (*Lijab)["mlab"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Fij)["lm"] * (*Lijab)["mkab"];
  (*LHijab)["klab"] += ( - 1.0  ) * (*Fab)["eb"] * (*Lijab)["klea"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Fab)["ea"] * (*Lijab)["kleb"];
  (*LHijab)["klab"] += ( - 0.5  ) * (*Vijkl)["klmn"] * (*Lijab)["nmab"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Viajb)["kenb"] * (*Lijab)["nlea"];
  (*LHijab)["klab"] += ( - 1.0  ) * (*Viajb)["kena"] * (*Lijab)["nleb"];
  (*LHijab)["klab"] += ( - 1.0  ) * (*Viajb)["lenb"] * (*Lijab)["nkea"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Viajb)["lena"] * (*Lijab)["nkeb"];
  (*LHijab)["klab"] += ( - 0.5  ) * (*Vabcd)["efab"] * (*Lijab)["klfe"];
  (*LHijab)["klab"] += ( + 0.5  ) * (*Tabij)["efop"] * (*Vijab)["klfb"] * (*Lijab)["poea"];
  (*LHijab)["klab"] += ( - 0.5  ) * (*Tabij)["efop"] * (*Vijab)["klfa"] * (*Lijab)["poeb"];
  (*LHijab)["klab"] += ( - 0.25  ) * (*Tabij)["efop"] * (*Vijab)["klef"] * (*Lijab)["poab"];
  (*LHijab)["klab"] += ( - 0.5  ) * (*Tabij)["efop"] * (*Vijab)["pkab"] * (*Lijab)["olfe"];
  (*LHijab)["klab"] += ( + 0.5  ) * (*Tabij)["efop"] * (*Vijab)["plab"] * (*Lijab)["okfe"];
  (*LHijab)["klab"] += ( - 1.0  ) * (*Tabij)["efop"] * (*Vijab)["pkfb"] * (*Lijab)["olea"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Tabij)["efop"] * (*Vijab)["pkfa"] * (*Lijab)["oleb"];
  (*LHijab)["klab"] += ( + 1.0  ) * (*Tabij)["efop"] * (*Vijab)["plfb"] * (*Lijab)["okea"];
  (*LHijab)["klab"] += ( - 1.0  ) * (*Tabij)["efop"] * (*Vijab)["plfa"] * (*Lijab)["okeb"];
  (*LHijab)["klab"] += ( + 0.5  ) * (*Tabij)["efop"] * (*Vijab)["pkef"] * (*Lijab)["olab"];
  (*LHijab)["klab"] += ( - 0.5  ) * (*Tabij)["efop"] * (*Vijab)["plef"] * (*Lijab)["okab"];
  (*LHijab)["klab"] += ( - 0.25  ) * (*Tabij)["efop"] * (*Vijab)["opab"] * (*Lijab)["klfe"];
  (*LHijab)["klab"] += ( - 0.5  ) * (*Tabij)["efop"] * (*Vijab)["opfb"] * (*Lijab)["klea"];
  (*LHijab)["klab"] += ( + 0.5  ) * (*Tabij)["efop"] * (*Vijab)["opfa"] * (*Lijab)["kleb"];

  return LH;
}

template <typename F>
void CcsdSimilarityTransformedHamiltonian<F>::buildIntermediates(
    bool intermediates, bool singleParticleOnly
  ) {

  withIntermediates = intermediates;

  if (! intermediates) {
    LOG(0, "CcsdEomDavid") << "Not building intermediates" << std::endl;
    return;
  }

  //[1]
  //Isaiah Shavitt, Rodney J. Bartlett. Many-Body Methods in Chemistry and
  //Physics: MBPT and Coupled-Cluster Theory. 2009
  //PAGE: 439

  //[2]
  //John F. Stanton, Rodney J. Bartlett. The equation of motion coupled‐cluster
  //method. A systematic biorthogonal approach to molecular excitation
  //energies, transition probabilities, and excited state properties. The
  //Journal of Chemical Physics 7029--7039  1993
  // TABLE 1

  LOG(0, "CcsdEomDavid") << "Building intermediates Wpqrs and Wpq"
                         << std::endl;
  auto Tau_abij(NEW(CTF::Tensor<>, *Tabij));
  (*Tau_abij)["abij"] += (*Tai)["ai"] * (*Tai)["bj"];
  (*Tau_abij)["abij"] += ( - 1.0 ) * (*Tai)["bi"] * (*Tai)["aj"];

  //This approach defines intermediates:
  //Wab Wia Wabcd Wabci Waibc Wiabj Wiajk Wij Wijka Wijkl

  int No(Fij->lens[0]);
  int Nv(Fab->lens[0]);
  int syms[] = {NS, NS};
  int ov[] = {No, Nv};
  CTF::Tensor<> InitFia(2, ov, syms, *Cc4s::world, "InitFia");

  Wia   = NEW(CTF::Tensor<>,  InitFia);
  Wab   = NEW(CTF::Tensor<>, *Fab);
  Wij   = NEW(CTF::Tensor<>, *Fij);
  Wabcd = NEW(CTF::Tensor<>, *Vabcd);
  Wabci = NEW(CTF::Tensor<>, *Vabci);
  Waibc = NEW(CTF::Tensor<>, *Vaibc);
  Wiabj = NEW(CTF::Tensor<>, *Viabj);
  Wiajk = NEW(CTF::Tensor<>, *Viajk);
  Wijka = NEW(CTF::Tensor<>, *Vijka);
  Wijkl = NEW(CTF::Tensor<>, *Vijkl);
  // Initialize intermediates to zero
  (*Wia)["do"] = 0.0;
  (*Wab)["do"] = 0.0;
  (*Wij)["do"] = 0.0;
  (*Wabcd)["blah"] = 0.0;
  (*Wabci)["blah"] = 0.0;
  (*Waibc)["blah"] = 0.0;
  (*Wiabj)["blah"] = 0.0;
  (*Wiajk)["blah"] = 0.0;
  (*Wijka)["blah"] = 0.0;
  (*Wijkl)["blah"] = 0.0;

  LOG(0, "CcsdEomDavid") << "Building Wia" << std::endl;
  //we need this one to construct the 2-body-amplitudes, not directly
  if (Fia) {
    (*Wia)["ia"] += (*Fia)["ia"];
  }
  if (!singleParticleOnly) {
    (*Wia)["ia"] = (*Vijab)["imae"] * (*Tai)["em"];
  }

  LOG(0, "CcsdEomDavid") << "Building Wab" << std::endl;
  //diagram (10.54)
  (*Wab)["ab"]  = (*Fab)["ab"];
  if (!singleParticleOnly) {
    (*Wab)["ab"] += (*Viabc)["mafb"] * (*Tai)["fm"];
    if (Fia) {
      (*Wab)["ab"] += ( -1.0) * (*Fia)["mb"] * (*Tai)["am"];
    }
    (*Wab)["ab"] += (- 0.5) * (*Vijab)["mnbe"] * (*Tau_abij)["aemn"];
  }

  LOG(0, "CcsdEomDavid") << "Building Wij" << std::endl;
  (*Wij)["ij"]  = (*Fij)["ij"];
  if (!singleParticleOnly) {
    (*Wij)["ij"] += (*Vijka)["imje"] * (*Tai)["em"];
    if (Fia) {
      (*Wij)["ij"] += (*Fia)["ie"] * (*Tai)["ej"];
    }
    (*Wij)["ij"] += (  0.5) * (*Vijab)["imef"] * (*Tau_abij)["efjm"];

    LOG(0, "CcsdEomDavid") << "Building Wijkl" << std::endl;
    //Taken directly from [2]
    (*Wijkl)["klij"]  = (*Vijkl)["klij"];
    //------------------------------------------------------------
    (*Wijkl)["klij"] +=           (*Tai)["ej"] * (*Vijka)["klie"];
    (*Wijkl)["klij"] += ( -1.0) * (*Tai)["ei"] * (*Vijka)["klje"];
    //------------------------------------------------------------
    (*Wijkl)["klij"] += ( 0.5 ) * (*Tau_abij)["efij"] * (*Vijab)["klef"];

    LOG(0, "CcsdEomDavid") << "Building Wabcd" << std::endl;
    (*Wabcd)["abcd"]  = (*Vabcd)["abcd"];
    //-----------------------------------------------------------
    (*Wabcd)["abcd"] += (-1.0) * (*Vaibc)["amcd"] * (*Tai)["bm"];
    // P(ab)
    (*Wabcd)["abcd"] += ( 1.0) * (*Vaibc)["bmcd"] * (*Tai)["am"];
    //-----------------------------------------------------------
    (*Wabcd)["abcd"] += ( 0.5) * (*Vijab)["mncd"] * (*Tau_abij)["abmn"];

    LOG(0, "CcsdEomDavid") << "Building Waibc" << std::endl;
    (*Waibc)["aibc"]  = (*Vaibc)["aibc"];
    (*Waibc)["aibc"] += ( -1.0) * (*Vijab)["mibc"] * (*Tai)["am"];

    LOG(0, "CcsdEomDavid") << "Building Wijka" << std::endl;
    //Taken directly from[2]
    (*Wijka)["jkia"]  = (*Vijka)["jkia"];
    (*Wijka)["jkia"] += (*Tai)["ei"] * (*Vijab)["jkea"];

    LOG(0, "CcsdEomDavid") << "Building Wiabj from Waijb" << std::endl;
    //[1] diagram (10.73)
    //This is not listed in the source book, however we can write it in terms
    //of Waijb since it should also have the simmetry of the Tabij amplitudes
    //and the Coulomb integrals Vpqrs
    //Taken directly from [2]
    (*Wiabj)["jabi"]  = (*Vaijb)["ajib"];
    (*Wiabj)["jabi"] += (*Vaibc)["ajeb"] * (*Tai)["ei"];
    (*Wiabj)["jabi"] += ( -1.0) * (*Vijka)["mjib"] * (*Tai)["am"];
    (*Wiabj)["jabi"] += ( -1.0) * (*Vijab)["mjeb"] * (*Tai)["ei"] * (*Tai)["am"];
    (*Wiabj)["jabi"] += ( -1.0) * (*Vijab)["mjeb"] * (*Tabij)["eaim"];

    bool wabciIntermediates(true);
    if (wabciIntermediates) {
      LOG(0, "CcsdEomDavid") << "Building Wabci from Wabcd and Wia" << std::endl;
      //--1
      (*Wabci)["abci"]  = (*Vabci)["abci"];
      //--3
      (*Wabci)["abci"] += ( -1.0) * (*Vaibj)["amci"] * (*Tai)["bm"];
      (*Wabci)["abci"] += ( +1.0) * (*Vaibj)["bmci"] * (*Tai)["am"];
      //--6
      (*Wabci)["abci"] += ( +1.0) * (*Vaibc)["amce"] * (*Tabij)["ebmi"];
      (*Wabci)["abci"] += ( -1.0) * (*Vaibc)["bmce"] * (*Tabij)["eami"];
      //--9
      (*Wabci)["abci"] += ( -1.0) * (*Vijab)["mnce"] * (*Tai)["am"] * (*Tabij)["ebni"];
      (*Wabci)["abci"] += ( +1.0) * (*Vijab)["mnce"] * (*Tai)["bm"] * (*Tabij)["eani"];
      //--8
      (*Wabci)["abci"] += ( -1.0) * (*Wia)["mc"] * (*Tabij)["abmi"];
      //--2-4-10-11
      (*Wabci)["abci"] += ( +1.0) * (*Tai)["ei"] * (*Wabcd)["abce"];
      //--7-5
      (*Wabci)["abci"] += (  0.5 ) * (*Vijak)["nmci"] * (*Tau_abij)["abnm"];
    } else {
      LOG(0, "CcsdEomDavid") << "Building Wabci" << std::endl;
      //--1
      (*Wabci)["abci"]  = (*Vabci)["abci"];
      //--2
      (*Wabci)["abci"] += (*Vabcd)["abce"] * (*Tai)["ei"];
      //--3
      (*Wabci)["abci"] += ( -1.0) * (*Vaibj)["amci"] * (*Tai)["bm"];
      (*Wabci)["abci"] += ( +1.0) * (*Vaibj)["bmci"] * (*Tai)["am"];
      //--4
      (*Wabci)["abci"] += ( -1.0) * (*Vaibc)["amce"] * (*Tai)["bm"] * (*Tai)["ei"];
      (*Wabci)["abci"] += ( +1.0) * (*Vaibc)["bmce"] * (*Tai)["am"] * (*Tai)["ei"];
      //--5
      (*Wabci)["abci"] += ( +1.0) * (*Vijak)["mnci"] * (*Tai)["am"] * (*Tai)["bn"];
      //--5.1 (non canonical)
      if (Fia) {
        (*Wabci)["abci"] += ( -1.0) * (*Fia)["mc"] * (*Tabij)["abmi"];
      }
      //--6
      (*Wabci)["abci"] +=          (*Vaibc)["amce"] * (*Tabij)["ebmi"];
      (*Wabci)["abci"] += (-1.0) * (*Vaibc)["bmce"] * (*Tabij)["eami"];
      //--7
      (*Wabci)["abci"] += (  0.5) * (*Vijak)["mnci"] * (*Tabij)["abmn"];
      //--8
      (*Wabci)["abci"] += ( -1.0) * (*Vijab)["mnec"] * (*Tai)["em"] * (*Tabij)["abni"];
      //--9
      (*Wabci)["abci"] += ( -1.0) * (*Vijab)["mnce"] * (*Tai)["am"] * (*Tabij)["ebni"];
      (*Wabci)["abci"] += ( +1.0) * (*Vijab)["mnce"] * (*Tai)["bm"] * (*Tabij)["eani"];
      //--10
      (*Wabci)["abci"] += (  0.5) * (*Vijab)["mnce"] * (*Tai)["ei"] * (*Tabij)["abmn"];
      //--11
      (*Wabci)["abci"] +=           (*Vijab)["mnce"] * (*Tai)["am"] * (*Tai)["bn"] * (*Tai)["ei"];
    }

    LOG(0, "CcsdEomDavid") << "Building Wiajk from Wia and Wijkl" << std::endl;
    //This is built upon the already existing amplitudes
    //[1] diagram (10.79)
    //Takend directly from [2]
    //--1
    (*Wiajk)["iajk"]  = (*Viajk)["iajk"];
    //--6
    (*Wiajk)["iajk"] +=            (*Vijka)["imje"] * (*Tabij)["aekm"];
    (*Wiajk)["iajk"] += ( -1.0 ) * (*Vijka)["imke"] * (*Tabij)["aejm"];
    //    original
    //    (*Wiajk)["iajk"] +=            (*Vijka)["imje"] * (*Tabij)["aekm"];
    //    (*Wiajk)["iajk"] += ( -1.0 ) * (*Vijka)["jmie"] * (*Tabij)["aekm"];
    //--7-5
    (*Wiajk)["iajk"] += (  0.5 ) * (*Viabc)["iaef"] * (*Tau_abij)["efjk"];
    //--8
    (*Wiajk)["iajk"] += ( -1.0) * (*Wia)["ie"] * (*Tabij)["aejk"];
    //    original: (Problem, The diagram actually says that it
    //    should be Teajk and not Taejk, so that 'a' stays in the second
    //    vertex, so we have to either change a<>e or put a minus)
    //    (*Wiajk)["iajk"] += (*Wia)["ie"] * (*Tabij)["aejk"];
    //--2-4-10-11
    (*Wiajk)["iajk"] += (-1.0) * (*Tai)["am"] * (*Wijkl)["imjk"];
    //    original: (minus)
    //    (*Wiajk)["iajk"] += (*Tai)["am"] * (*Wijkl)["imjk"];
    //--3
    (*Wiajk)["iajk"] += ( +1.0 ) * (*Tai)["ek"] * (*Viajb)["iaje"];
    (*Wiajk)["iajk"] += ( -1.0 ) * (*Tai)["ej"] * (*Viajb)["iake"];
    //     original:
    //     (*Wiajk)["iajk"] += ( -1.0 ) * (*Tai)["ej"] * (*Viabj)["iaek"];
    //     (*Wiajk)["iajk"] += ( +1.0 ) * (*Tai)["ei"] * (*Viabj)["jaek"];
    //--9
    (*Wiajk)["iajk"] +=
      ( -1.0 ) * (*Tai)["ej"] * (*Tabij)["afmk"] * (*Vijab)["imef"];
    (*Wiajk)["iajk"] +=
      ( +1.0 ) * (*Tai)["ek"] * (*Tabij)["afmj"] * (*Vijab)["imef"];
    //     original: Again it does not make any sense to do Pij, and the minus
    //     (*Wiajk)["iajk"] +=
    //       ( +1.0 ) * (*Tai)["ej"] * (*Tabij)["afmk"] * (*Vijab)["imef"];
    //     (*Wiajk)["iajk"] +=
    //       ( -1.0 ) * (*Tai)["ei"] * (*Tabij)["afmk"] * (*Vijab)["jmef"];
  }
}

template <typename F>
FockVector<F> CcsdSimilarityTransformedHamiltonian<F>::rightApply(
  FockVector<F> &R
) {
  return withIntermediates ? rightApplyIntermediates(R) : rightApplyHirata(R);
}

template <typename F>
FockVector<F> CcsdSimilarityTransformedHamiltonian<F>::rightApplyHirata(
  FockVector<F> &R
) {
  FockVector<F> HR(R);
  // get pointers to the component tensors
  PTR(CTF::Tensor<F>) Rai( R.get(0) );
  PTR(CTF::Tensor<F>) Rabij( R.get(1) );
  PTR(CTF::Tensor<F>) HRai( HR.get(0) );
  PTR(CTF::Tensor<F>) HRabij( HR.get(1) );

  checkAntisymmetry(*Rabij);

  // Contruct HR (one body part)
  // TODO: why "bi" not "ai"?
  (*HRai)["bi"]  = 0.0;

  // WIJ =====================================================================
  (*HRai)["bi"] += ( - 1.0 ) * (*Fij)["ki"] * (*Rai)["bk"];
  (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["cl"] * (*Vijka)["lmic"] * (*Rai)["bm"];
  (*HRai)["bi"] += ( - 0.5 ) * (*Tabij)["cdmi"] * (*Vijab)["mncd"] * (*Rai)["bn"];
  (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["ci"] * (*Tai)["dm"] * (*Vijab)["mncd"] * (*Rai)["bn"];

  // WAB =====================================================================
  (*HRai)["bi"] += ( + 1.0 ) * (*Fab)["bc"] * (*Rai)["ci"];
  (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["cl"] * (*Viabc)["lbce"] * (*Rai)["ei"];
  (*HRai)["bi"] += ( - 0.5 ) * (*Tabij)["cblm"] * (*Vijab)["lmcf"] * (*Rai)["fi"];
  (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["bk"] * (*Tai)["dm"] * (*Vijab)["kmdf"] * (*Rai)["fi"];

  // WIABJ ===================================================================
  (*HRai)["bi"] += ( - 1.0 ) * (*Viajb)["kbid"] * (*Rai)["dk"];
  (*HRai)["bi"] += ( + 1.0 ) * (*Tabij)["cbli"] * (*Vijab)["lmcf"] * (*Rai)["fm"];
  (*HRai)["bi"] += ( - 1.0  ) * (*Tai)["bk"] * (*Vijka)["klie"] * (*Rai)["el"];
  (*HRai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Viabc)["lbce"] * (*Rai)["el"];
  (*HRai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["bl"] * (*Vijab)["lmcf"] * (*Rai)["fm"];

  // WIA =====================================================================
  (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["cl"] * (*Vijab)["lmcf"] * (*Rabij)["fbmi"];

  // WIJKA ===================================================================
  (*HRai)["bi"] += ( + 0.5 ) * (*Vijka)["klie"] * (*Rabij)["ebkl"];
  (*HRai)["bi"] += ( + 0.5  ) * (*Tai)["ci"] * (*Vijab)["lmcf"] * (*Rabij)["fblm"];

  // WIABC ===================================================================
  (*HRai)["bi"] += ( + 0.5 ) * (*Viabc)["kbde"] * (*Rabij)["deki"];
  (*HRai)["bi"] += ( + 0.5  ) * (*Tai)["bk"] * (*Vijab)["klef"] * (*Rabij)["efli"];


  //(*HRai)["ai"]  = 0.0; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  (*HRabij)["cdij"]  = 0.0;

  // Contruct HR (two body part)

  // WABCD ===================================================================
  (*HRabij)["cdij"] += ( + 0.5 ) * (*Vabcd)["cdef"] * (*Rabij)["efij"];
  (*HRabij)["cdij"] +=
    ( - 0.5  ) * (*Tai)["cm"] * (*Viabc)["mdfg"] * (*Rabij)["fgij"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tai)["dm"] * (*Viabc)["mcfg"] * (*Rabij)["fgij"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijab)["mngh"] * (*Rabij)["ghij"];
  (*HRabij)["cdij"] +=
    ( + 0.25) * (*Tabij)["cdmn"] * (*Vijab)["mngh"] * (*Rabij)["ghij"];

  // WIJKL ===================================================================
  (*HRabij)["cdij"] += ( + 0.5 ) * (*Vijkl)["mnij"] * (*Rabij)["cdmn"];
  (*HRabij)["cdij"] +=
    ( + 0.25) * (*Tabij)["efij"] * (*Vijab)["opef"] * (*Rabij)["cdop"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tai)["ej"] * (*Vijka)["noie"] * (*Rabij)["cdno"];
  (*HRabij)["cdij"] +=
    ( - 0.5  ) * (*Tai)["ei"] * (*Vijka)["noje"] * (*Rabij)["cdno"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Vijab)["opef"] * (*Rabij)["cdop"];

  // WAB   ===================================================================
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Fab)["de"] * (*Rabij)["ecij"];
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Fab)["ce"] * (*Rabij)["edij"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["en"] * (*Viabc)["ndeg"] * (*Rabij)["gcij"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["en"] * (*Viabc)["nceg"] * (*Rabij)["gdij"];
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["edno"] * (*Vijab)["noeh"] * (*Rabij)["hcij"];
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["ecno"] * (*Vijab)["noeh"] * (*Rabij)["hdij"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["dm"] * (*Tai)["fo"] * (*Vijab)["mofh"] * (*Rabij)["hcij"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["fo"] * (*Vijab)["mofh"] * (*Rabij)["hdij"];

  // WIJ   ===================================================================
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Fij)["mi"] * (*Rabij)["cdmj"];
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Fij)["mj"] * (*Rabij)["cdmi"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["en"] * (*Vijka)["noie"] * (*Rabij)["cdoj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["en"] * (*Vijka)["noje"] * (*Rabij)["cdoi"];
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["efoi"] * (*Vijab)["opef"] * (*Rabij)["cdpj"];
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["efoj"] * (*Vijab)["opef"] * (*Rabij)["cdpi"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fo"] * (*Vijab)["opef"] * (*Rabij)["cdpj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["fo"] * (*Vijab)["opef"] * (*Rabij)["cdpi"];

  // WIABJ ===================================================================
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Viajb)["mdif"] * (*Rabij)["fcmj"];
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Viajb)["mcif"] * (*Rabij)["fdmj"];
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Viajb)["mdjf"] * (*Rabij)["fcmi"];
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Viajb)["mcjf"] * (*Rabij)["fdmi"];
  //--
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["dm"] * (*Vijka)["mnig"] * (*Rabij)["gcnj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["cm"] * (*Vijka)["mnig"] * (*Rabij)["gdnj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["dm"] * (*Vijka)["mnjg"] * (*Rabij)["gcni"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["cm"] * (*Vijka)["mnjg"] * (*Rabij)["gdni"];
  //--
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Viabc)["ndeg"] * (*Rabij)["gcnj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Viabc)["nceg"] * (*Rabij)["gdnj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Viabc)["ndeg"] * (*Rabij)["gcni"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ej"] * (*Viabc)["nceg"] * (*Rabij)["gdni"];
  //--
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Vijab)["noeh"] * (*Rabij)["hcoj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Vijab)["noeh"] * (*Rabij)["hdoj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Vijab)["noeh"] * (*Rabij)["hcoi"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Vijab)["noeh"] * (*Rabij)["hdoi"];
  //--
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["edni"] * (*Vijab)["noeh"] * (*Rabij)["hcoj"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ecni"] * (*Vijab)["noeh"] * (*Rabij)["hdoj"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ednj"] * (*Vijab)["noeh"] * (*Rabij)["hcoi"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["ecnj"] * (*Vijab)["noeh"] * (*Rabij)["hdoi"];

  //THREE_BODY_ONE ===========================================================
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ecij"] * (*Viabc)["ndeg"] * (*Rai)["gn"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["edij"] * (*Viabc)["nceg"] * (*Rai)["gn"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["ecij"] * (*Tai)["dn"] * (*Vijab)["noeh"] * (*Rai)["ho"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["edij"] * (*Tai)["cn"] * (*Vijab)["noeh"] * (*Rai)["ho"];

  //THREE_BODY_TWO ===========================================================
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["edij"] * (*Vijab)["noeh"] * (*Rabij)["hcno"];
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["ecij"] * (*Vijab)["noeh"] * (*Rabij)["hdno"];

  //THREE_BODY_THREE =========================================================
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["cdmj"] * (*Vijka)["mnig"] * (*Rai)["gn"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["cdmi"] * (*Vijka)["mnjg"] * (*Rai)["gn"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["cdmj"] * (*Tai)["fi"] * (*Vijab)["mofh"] * (*Rai)["ho"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["cdmi"] * (*Tai)["fj"] * (*Vijab)["mofh"] * (*Rai)["ho"];

  //THREE_BODY_FOUR ==========================================================
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["cdmi"] * (*Vijab)["mngh"] * (*Rabij)["ghnj"];
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["cdmj"] * (*Vijab)["mngh"] * (*Rabij)["ghni"];

  // WIAJK ===================================================================
  //--1
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Viajk)["mdij"] * (*Rai)["cm"];
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Viajk)["mcij"] * (*Rai)["dm"];
  //--2
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["dm"] * (*Vijkl)["mnij"] * (*Rai)["cn"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["cm"] * (*Vijkl)["mnij"] * (*Rai)["dn"];
  //--3
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Viajb)["ndie"] * (*Rai)["cn"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ej"] * (*Viajb)["ncie"] * (*Rai)["dn"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Viajb)["ndje"] * (*Rai)["cn"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Viajb)["ncje"] * (*Rai)["dn"];
  //--4
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Vijka)["noie"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Vijka)["noie"] * (*Rai)["do"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Vijka)["noje"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Vijka)["noje"] * (*Rai)["do"];
  //--5
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Viabc)["odef"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Viabc)["ocef"] * (*Rai)["do"];
  //--6
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ednj"] * (*Vijka)["noie"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["ecnj"] * (*Vijka)["noie"] * (*Rai)["do"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["edni"] * (*Vijka)["noje"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ecni"] * (*Vijka)["noje"] * (*Rai)["do"];
  //--7
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["efij"] * (*Viabc)["odef"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["efij"] * (*Viabc)["ocef"] * (*Rai)["do"];
  //--8
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["edij"] * (*Tai)["fo"] * (*Vijab)["opef"] * (*Rai)["cp"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["ecij"] * (*Tai)["fo"] * (*Vijab)["opef"] * (*Rai)["dp"];
  //--9
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["ednj"] * (*Tai)["gi"] * (*Vijab)["npeg"] * (*Rai)["cp"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["gi"] * (*Vijab)["npeg"] * (*Rai)["dp"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["edni"] * (*Tai)["gj"] * (*Vijab)["npeg"] * (*Rai)["cp"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["ecni"] * (*Tai)["gj"] * (*Vijab)["npeg"] * (*Rai)["dp"];
  //--10
  (*HRabij)["cdij"] +=
    ( - 0.5  ) * (*Tabij)["efij"] * (*Tai)["do"] * (*Vijab)["opef"] * (*Rai)["cp"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tabij)["efij"] * (*Tai)["co"] * (*Vijab)["opef"] * (*Rai)["dp"];
  //--11
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["do"] * (*Vijab)["opef"] * (*Rai)["cp"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Vijab)["opef"] * (*Rai)["dp"];

  // WABCI ===================================================================
  //--1
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Vabic)["cdie"] * (*Rai)["ej"];
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Vabic)["cdje"] * (*Rai)["ei"];
  //--2
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Vabcd)["cdef"] * (*Rai)["fj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Vabcd)["cdef"] * (*Rai)["fi"];
  //--3
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["cm"] * (*Viajb)["mdif"] * (*Rai)["fj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["dm"] * (*Viajb)["mcif"] * (*Rai)["fj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["cm"] * (*Viajb)["mdjf"] * (*Rai)["fi"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["dm"] * (*Viajb)["mcjf"] * (*Rai)["fi"];
  //--4
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Viabc)["ndeg"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Viabc)["nceg"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Viabc)["ndeg"] * (*Rai)["gi"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Viabc)["nceg"] * (*Rai)["gi"];
  //--5
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijka)["mnig"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijka)["mnjg"] * (*Rai)["gi"];
  //--6
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ecni"] * (*Viabc)["ndeg"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["edni"] * (*Viabc)["nceg"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["ecnj"] * (*Viabc)["ndeg"] * (*Rai)["gi"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ednj"] * (*Viabc)["nceg"] * (*Rai)["gi"];
  //--7
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["cdmn"] * (*Vijka)["mnig"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["cdmn"] * (*Vijka)["mnjg"] * (*Rai)["gi"];
  //--8
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["cdmi"] * (*Tai)["fo"] * (*Vijab)["mofh"] * (*Rai)["hj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["cdmj"] * (*Tai)["fo"] * (*Vijab)["mofh"] * (*Rai)["hi"];
  //--9
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["ecni"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["edni"] * (*Tai)["co"] * (*Vijab)["noeh"] * (*Rai)["hj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hi"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["ednj"] * (*Tai)["co"] * (*Vijab)["noeh"] * (*Rai)["hi"];
  //--10
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gi"] * (*Vijab)["mngh"] * (*Rai)["hj"];
  (*HRabij)["cdij"] +=
    ( - 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gj"] * (*Vijab)["mngh"] * (*Rai)["hi"];
  //--11
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hj"];
  (*HRabij)["cdij"] +=
   ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hi"];

  // NON CANONICAL ORBITALS ==================================================

  if ( Fia ) {

    (*HRai)["bi"] +=  ( + 1.0  ) * (*Fia)["kd"] * (*Rabij)["dbki"];
    (*HRai)["bi"] +=  ( - 1.0  ) * (*Fia)["kd"] * (*Tai)["di"] * (*Rai)["bk"];
    (*HRai)["bi"] +=  ( - 1.0  ) * (*Fia)["kd"] * (*Tai)["bk"] * (*Rai)["di"];

    (*HRabij)["cdij"] += ( - 1.0  ) * (*Fia)["mf"] * (*Tai)["fi"] * (*Rabij)["cdmj"];
    (*HRabij)["cdij"] += ( + 1.0  ) * (*Fia)["mf"] * (*Tai)["fj"] * (*Rabij)["cdmi"];

    (*HRabij)["cdij"] += ( + 1.0  ) * (*Fia)["mf"] * (*Tai)["dm"] * (*Rabij)["fcij"];
    (*HRabij)["cdij"] += ( - 1.0  ) * (*Fia)["mf"] * (*Tai)["cm"] * (*Rabij)["fdij"];

    (*HRabij)["cdij"] += ( - 1.0  ) * (*Fia)["mf"] * (*Tabij)["fdij"] * (*Rai)["cm"];
    (*HRabij)["cdij"] += ( + 1.0  ) * (*Fia)["mf"] * (*Tabij)["fcij"] * (*Rai)["dm"];

    (*HRabij)["cdij"] += ( + 1.0  ) * (*Fia)["mf"] * (*Tabij)["cdmi"] * (*Rai)["fj"];
    (*HRabij)["cdij"] += ( - 1.0  ) * (*Fia)["mf"] * (*Tabij)["cdmj"] * (*Rai)["fi"];

  }

  return HR;

}

template <typename F>
FockVector<F> CcsdSimilarityTransformedHamiltonian<F>::rightApplyIntermediates(
  FockVector<F> &R
) {
  FockVector<F> HR(R);
  // get pointers to the component tensors
  PTR(CTF::Tensor<F>) Rai( R.get(0) );
  PTR(CTF::Tensor<F>) Rabij( R.get(1) );
  PTR(CTF::Tensor<F>) HRai( HR.get(0) );
  PTR(CTF::Tensor<F>) HRabij( HR.get(1) );

  checkAntisymmetry(*Rabij);

  (*HRai)["ai"]  = 0.0;
  (*HRai)["ai"] += (- 1.0) * (*Wij)["li"] * (*Rai)["al"];
  (*HRai)["ai"] += (*Wab)["ad"] * (*Rai)["di"];
  (*HRai)["ai"] += (*Wiabj)["ladi"] * (*Rai)["dl"];

  (*HRai)["ai"] += (*Wia)["ld"] * (*Rabij)["adil"];

  (*HRai)["ai"] += ( - 0.5 ) * (*Wijka)["lmid"] * (*Rabij)["adlm"];
  (*HRai)["ai"] += (   0.5 ) * (*Waibc)["alde"] * (*Rabij)["deil"];

  //(*HRai)["ai"]  = 0.0; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // 2 body part
  (*HRabij)["abij"]  = 0.0;

  // WABCD ===================================================================
  (*HRabij)["abij"] += (  0.5) * (*Wabcd)["abde"] * (*Rabij)["deij"];

  // WIJKL ===================================================================
  (*HRabij)["abij"] += (  0.5) * (*Wijkl)["lmij"] * (*Rabij)["ablm"];

  // WAB   ===================================================================
  (*HRabij)["abij"] += ( +1.0) * (*Wab)["bd"] * (*Rabij)["adij"];
  //P(ab)
  (*HRabij)["abij"] += ( -1.0) * (*Wab)["ad"] * (*Rabij)["bdij"];

  // WIJ   ===================================================================
  (*HRabij)["abij"] += ( -1.0) * (*Wij)["lj"] * (*Rabij)["abil"];
  //P(ij)
  (*HRabij)["abij"] +=           (*Wij)["li"] * (*Rabij)["abjl"];

  // WIABJ ===================================================================
  (*HRabij)["abij"] +=            (*Wiabj)["lbdj"] * (*Rabij)["adil"];
  //-P(ij)
  (*HRabij)["abij"] +=  ( -1.0) * (*Wiabj)["lbdi"] * (*Rabij)["adjl"];
  //-P(ab)
  (*HRabij)["abij"] +=  ( -1.0) * (*Wiabj)["ladj"] * (*Rabij)["bdil"];
  //P(ij)P(ab)
  (*HRabij)["abij"] +=            (*Wiabj)["ladi"] * (*Rabij)["bdjl"];

  //THREE_BODY_ONE ===========================================================
  (*HRabij)["abij"] +=
              (*Rai)["em"] * (*Waibc)["bmfe"] * (*Tabij)["afij"];
  // P(ab)
  (*HRabij)["abij"] +=
    ( -1.0) * (*Rai)["em"] * (*Waibc)["amfe"] * (*Tabij)["bfij"];

  //THREE_BODY_TWO ===========================================================
  (*HRabij)["abij"] +=
    ( -0.5) * (*Rabij)["eamn"] * (*Vijab)["nmfe"] * (*Tabij)["fbij"];
  // P(ab)
  (*HRabij)["abij"] +=
    ( +0.5) * (*Rabij)["ebmn"] * (*Vijab)["nmfe"] * (*Tabij)["faij"];

  //THREE_BODY_THREE =========================================================
  (*HRabij)["abij"] +=
    ( -1.0) * (*Rai)["em"] * (*Wijka)["nmje"] * (*Tabij)["abin"];
  // P(ij)
  (*HRabij)["abij"] +=
    ( +1.0) * (*Rai)["em"] * (*Wijka)["nmie"] * (*Tabij)["abjn"];

  //THREE_BODY_FOUR ==========================================================
  (*HRabij)["abij"] +=
    ( +0.5) * (*Rabij)["feim"] * (*Vijab)["nmfe"] * (*Tabij)["abjn"];
  // P(ij)
  (*HRabij)["abij"] +=
    ( -0.5) * (*Rabij)["fejm"] * (*Vijab)["nmfe"] * (*Tabij)["abin"];

  // WIAJK ===================================================================
  (*HRabij)["abij"] += (- 1.0 ) * (*Wiajk)["lbij"] * (*Rai)["al"];
  //P(ab)
  (*HRabij)["abij"] += (+ 1.0 ) * (*Wiajk)["laij"] * (*Rai)["bl"];

  // WABCI ===================================================================
  (*HRabij)["abij"] +=          (*Wabci)["abej"] * (*Rai)["ei"];
  //P(ij)
  (*HRabij)["abij"] += (-1.0) * (*Wabci)["abei"] * (*Rai)["ej"];

  return HR;
}

// instantiate class
template
class CcsdSimilarityTransformedHamiltonian<double>;


template <typename F>
CcsdPreConditioner<F>::CcsdPreConditioner(
  CTF::Tensor<F> &Tai,
  CTF::Tensor<F> &Tabij,
  CTF::Tensor<F> &Fij,
  CTF::Tensor<F> &Fab,
  CTF::Tensor<F> &Vabcd,
  CTF::Tensor<F> &Viajb,
  CTF::Tensor<F> &Vijab,
  CTF::Tensor<F> &Vijkl,
  Algorithm *algorithm_
): diagonalH(
    std::vector<PTR(CTF::Tensor<double>)>(
      {NEW(CTF::Tensor<double>, Tai), NEW(CTF::Tensor<double>, Tabij)}
    ),
    std::vector<std::string>({"ai", "abij"})
  ), algorithm(algorithm_)
  {
  // pointers to singles and doubles tensors of diagonal part
  auto Dai( diagonalH.get(0) );
  auto Dabij( diagonalH.get(1) );

  // TODO: Maybe inster the Tai part to the diagonal

  // calculate diagonal elements of H
  (*Dai)["bi"] =  ( - 1.0 ) * Fij["ii"];
  (*Dai)["bi"] += ( + 1.0 ) * Fab["bb"];
  (*Dai)["bi"] += ( - 1.0 ) * Viajb["ibib"];
  (*Dai)["bi"] += ( + 1.0 ) * Tabij["cbli"] * Vijab["licb"];
  (*Dai)["bi"] += ( - 0.5 ) * Tabij["cdmi"] * Vijab["micd"];
  (*Dai)["bi"] += ( - 0.5 ) * Tabij["cblm"] * Vijab["lmcb"];

  (*Dabij)["cdij"] =  ( - 1.0 ) * Fij["ii"];
  (*Dabij)["cdii"] += ( + 1.0 ) * Fij["ii"];
  (*Dabij)["ccij"] += ( - 1.0 ) * Fab["cc"];
  (*Dabij)["cdij"] += ( + 1.0 ) * Fab["cc"];

  (*Dabij)["cdij"] += ( + 0.5 ) * Vijkl["ijij"];
  (*Dabij)["ccij"] += ( + 1.0 ) * Viajb["icic"];
  (*Dabij)["cdij"] += ( - 1.0 ) * Viajb["icic"];
  (*Dabij)["ccii"] += ( - 1.0 ) * Viajb["icic"];
  (*Dabij)["cdii"] += ( + 1.0 ) * Viajb["icic"];
  (*Dabij)["cdij"] += ( + 0.5 ) * Vabcd["cdcd"];
  (*Dabij)["ccij"] += ( + 0.5 ) * Tabij["ecij"] * Vijab["ijec"];
  (*Dabij)["cdij"] += ( - 0.5 ) * Tabij["ecij"] * Vijab["ijec"];
  (*Dabij)["cdij"] += ( + 0.25) * Tabij["efij"] * Vijab["ijef"];
  (*Dabij)["cdij"] += ( - 0.5 ) * Tabij["cdmi"] * Vijab["micd"];
  (*Dabij)["cdii"] += ( + 0.5 ) * Tabij["cdmi"] * Vijab["micd"];
  (*Dabij)["ccij"] += ( - 1.0 ) * Tabij["ecni"] * Vijab["niec"];
  (*Dabij)["cdij"] += ( + 1.0 ) * Tabij["ecni"] * Vijab["niec"];
  (*Dabij)["ccii"] += ( + 1.0 ) * Tabij["ecni"] * Vijab["niec"];
  (*Dabij)["cdii"] += ( - 1.0 ) * Tabij["ecni"] * Vijab["niec"];
  (*Dabij)["cdij"] += ( - 0.5 ) * Tabij["efoi"] * Vijab["oief"];
  (*Dabij)["cdii"] += ( + 0.5 ) * Tabij["efoi"] * Vijab["oief"];
  (*Dabij)["cdij"] += ( + 0.25) * Tabij["cdmn"] * Vijab["mncd"];
  (*Dabij)["ccij"] += ( + 0.5 ) * Tabij["ecno"] * Vijab["noec"];
  (*Dabij)["cdij"] += ( - 0.5 ) * Tabij["ecno"] * Vijab["noec"];

}

template <typename F>
class EomDiagonalValueComparator;

/**
 * \brief Comparator that should filter out zero values of the diagonal
 * matrix.
 * Zero values are treated as infinite so that they get appended to the
 * end of the list.
 */
template <>
class EomDiagonalValueComparator<double> {
public:
  bool operator ()(
    const std::pair<int, double> &a,
    const std::pair<int, double> &b
  ) {
    double A(
      std::abs(a.second) < 1E-13 ?
        std::numeric_limits<double>::infinity() : a.second
    );
    double B(
      std::abs(b.second) < 1E-13 ?
        std::numeric_limits<double>::infinity() : b.second
    );
    double diff(B-A);
    // maintain magnitude finite!
    double magnitude(std::abs(a.second)+std::abs(b.second));
    if (std::real(diff) > +1E-13*magnitude) return true;
    if (std::real(diff) < -1E-13*magnitude) return false;
    return a.first < b.first;
  }
};

template <typename F>
std::vector<FockVector<F>>
CcsdPreConditioner<F>::getInitialBasis(
  const int eigenVectorsCount
) {
  LOG(0, "CcsdPreConditioner") << "Getting initial basis " << std::endl;
  int random(algorithm->getIntegerArgument("preconditionerRandom", 0));
  DefaultRandomEngine randomEngine;
  std::normal_distribution<double> normalDistribution(
    0.0, algorithm->getRealArgument("preconditionerRandomSigma", 0.1)
  );
  if (random == 1) {
    LOG(0, "CcsdPreConditioner") << "Randomizing initial guess" << std::endl;
  }
  // find K=eigenVectorsCount lowest diagonal elements at each processor
  std::vector<std::pair<size_t, F>> localElements( diagonalH.readLocal() );
  std::sort(
    localElements.begin(), localElements.end(),
    EomDiagonalValueComparator<double>()
  );
  int localElementsSize( localElements.size() );

  // gather all K elements of all processors at root
  //   convert into homogeneous arrays for MPI gather
  const int trialEigenVectorsCount(10*eigenVectorsCount);
  std::vector<size_t> localLowestElementIndices(trialEigenVectorsCount);
  std::vector<F> localLowestElementValues(trialEigenVectorsCount);
  for (
    size_t i(0);
    i < std::min(localElementsSize, trialEigenVectorsCount);
    ++i
  ) {
    localLowestElementIndices[i] = localElements[i].first;
    localLowestElementValues[i] = localElements[i].second;
  }
  MpiCommunicator communicator(*Cc4s::world);
  std::vector<size_t> lowestElementIndices;
  std::vector<F> lowestElementValues;
  communicator.gather(localLowestElementIndices, lowestElementIndices);
  communicator.gather(localLowestElementValues, lowestElementValues);
  //   convert back into (index,value) pairs for sorting
  std::vector<std::pair<size_t, F>> lowestElements(lowestElementValues.size());
  for (int i(0); i < lowestElementValues.size(); ++i) {
    lowestElements[i].first = lowestElementIndices[i];
    lowestElements[i].second = lowestElementValues[i];
  }

  // find globally lowest K diagonal elements among the gathered elements
  std::sort(
    lowestElements.begin(), lowestElements.end(),
    EomDiagonalValueComparator<double>()
  );
  // at rank==0 (root) lowestElements contains K*Np entries
  // rank > 0 has an empty list

  // create basis vectors for each lowest element
  std::vector<V> basis;

  int currentEigenVectorCount(0);
  int b(0);
  int zeroVectorCount(0);
  while (currentEigenVectorCount < eigenVectorsCount) {
    V basisElement(diagonalH);
    basisElement *= 0.0;
    std::vector<std::pair<size_t,F>> elements;
    if (communicator.getRank() == 0) {
      elements.push_back(
        std::make_pair(lowestElements[b].first, 1.0)
      );
    }
    basisElement.write(elements);
    if (random == 1) {
      auto Rai(*basisElement.get(0));
      auto Rabij(*basisElement.get(1));
      setRandomTensor(Rai, normalDistribution, randomEngine);
      setRandomTensor(Rabij, normalDistribution, randomEngine);
      (*basisElement.get(0))["ai"] += Rai["ai"];
      (*basisElement.get(1))["abij"] += Rabij["abij"];
    }
    // (101, -70), (32, -55), ...
    // b1: 0... 1 (at global position 101) 0 ...
    // b2: 0... 1 (at global position 32) 0 ...i

    // Filter out unphysical components from the basisElement
    (*basisElement.get(1))["abii"]=0.0;
    (*basisElement.get(1))["aaij"]=0.0;
    (*basisElement.get(1))["aaii"]=0.0;

    // Antisymmetrize the new basis element
    (*basisElement.get(1))["abij"] -= (*basisElement.get(1))["abji"];
    (*basisElement.get(1))["abij"] -= (*basisElement.get(1))["baij"];

    // Grams-schmidt it with the other elements of the basis
    for (unsigned int j(0); j < basis.size(); ++j) {
      basisElement -= basis[j] * basisElement.dot(basis[j]);
    }

    // Normalize basisElement
    F basisElementNorm(std::sqrt(basisElement.dot(basisElement)));

    // Check if basisElementNorm is zero
    if ( std::abs(basisElementNorm) < 1e-10 ) {
      zeroVectorCount++;
      b++;
      continue;
    }

    basisElement = 1 / std::sqrt(basisElement.dot(basisElement))*basisElement;
    basisElementNorm = std::sqrt(basisElement.dot(basisElement));

    b++;

    if ( std::abs(basisElementNorm - F(1)) > 1e-10 * F(1)) continue;

    currentEigenVectorCount++;

    // If it got here, basisElement is a valid vector
    basis.push_back(basisElement);

  }

  // Now make sure that you are giving an antisymmetric basis
  // and also that it is grammschmited afterwards
  //LOG(1,"CcsdPreConditioner") << "Antisymmetrize basis" << std::endl;
  //for (unsigned int j(0); j < basis.size(); ++j) {
    //(*basis[j].get(1))["abij"] -= (*basis[j].get(1))["abji"];
    //(*basis[j].get(1))["abij"] -= (*basis[j].get(1))["baij"];
  //}

  //LOG(1,"CcsdPreConditioner") <<
      //"Performing Gramm Schmidt in the initial basis"
  //<< std::endl;
  //for (unsigned int b(0); b < basis.size(); ++b) {
    //V newVector(basis[b]);
    //for (unsigned int j(0); j < b; ++j) {
      //newVector -= basis[j] * basis[j].dot(basis[b]);
    //}
    //// normalize
    //basis[b] = 1 / std::sqrt(newVector.dot(newVector)) * newVector;
  //}

  return basis;
}

template <typename F>
FockVector<F>
CcsdPreConditioner<F>::getCorrection(
  const complex lambda, FockVector<F> &residuum
) {
  FockVector<F> w(diagonalH);

  // Define a singleton helping class for the diagonal correction
  class DiagonalCorrection {
    public:
      DiagonalCorrection(const double lambda_): lambda(lambda_) {
      }
      F operator ()(const F residuumElement, const F diagonalElement) {
        return std::abs(lambda - diagonalElement) < 1E-4 ?
          0.0 : residuumElement / (lambda - diagonalElement);
      }
    protected:
      double lambda;
  } diagonalCorrection(std::real(lambda));

  FockVector<F> correction(diagonalH);
  // compute ((lambda * id - Diag(diagonal))^-1) . residuum
  for (unsigned int c(0); c < w.getComponentsCount(); ++c) {
    const char *indices( correction.componentIndices[c].c_str() );
    (*correction.get(c)).contract(
      1.0,
      *residuum.get(c),indices,
      *diagonalH.get(c),indices,
      0.0,indices,
      CTF::Bivar_Function<F>(diagonalCorrection)
    );
  }
  // Filter out unphysical components from the correction
  (*correction.get(1))["abii"]=0.0;
  (*correction.get(1))["aaij"]=0.0;
  (*correction.get(1))["aaii"]=0.0;

  // Antisymmetrize the correction
  (*correction.get(1))["abij"] -= (*correction.get(1))["abji"];
  (*correction.get(1))["abij"] -= (*correction.get(1))["baij"];
  (*correction.get(1))["abij"] = 0.25 * (*correction.get(1))["abij"];

  return correction;
}


// instantiate class
template
class CcsdPreConditioner<double>;

