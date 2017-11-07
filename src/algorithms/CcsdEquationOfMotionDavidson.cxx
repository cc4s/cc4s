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
#include <util/Exception.hpp>
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
  int kineticLensVirtual[] = {Nv, Nv};
  int kineticSyms[] = {NS, NS};
  CTF::Tensor<> Fab(2, kineticLensVirtual, kineticSyms, *Cc4s::world, "Fab");
  int kineticLensOccupied[] = {No, No};
  CTF::Tensor<> Fij(2, kineticLensOccupied, kineticSyms, *Cc4s::world, "Fij");
  Fab["aa"] = (*epsa)["a"];
  Fij["ii"] = (*epsi)["i"];

  // Get the Uccsd amplitudes
  CTF::Tensor<> Tai(
      getTensorArgument<double, CTF::Tensor<> >("SinglesAmplitudes"));
  CTF::Tensor<> Tabij(
      getTensorArgument<double, CTF::Tensor<> >("DoublesAmplitudes"));

  CcsdSimilarityTransformedHamiltonian<double> H(
    &Tai, &Tabij, &Fij, &Fab,
    Vabcd, Viajb, Vijab, Vijkl, Vijka, Viabc, Viajk, Vabic,
    Vaibc, Vaibj, Viabj, Vijak, Vaijb, Vabci
  );
  CcsdPreConditioner<double> P(
    Tai, Tabij, Fij, Fab, *Vabcd, *Viajb, *Vijab, *Vijkl
  );
  std::vector<FockVector<double>> basis( P.getInitialBasis(4) );
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
  LOG(0, "CcsdEomDavid") << "Computing " << eigenStates << " eigen states"
                              << std::endl;
  EigenSystemDavidson<FockVector<double>> eigenSystem(H, eigenStates, P, 1E-4, 8*16);

  std::vector<complex> eigenValues(eigenSystem.getEigenValues());
  for (auto &ev: eigenValues) {
    LOG(0, "CcsdEomDavid") << "Eigenvalue=" << ev << std::endl;
  }
}

// template method implementation
template <typename F>
void CcsdEquationOfMotionDavidson::getCanonicalPerturbationBasis(
  CTF::Tensor<F> &Tai, CTF::Tensor<F> &Tabij, size_t i
) {
  std::vector<std::pair<size_t, F>> elements;
  if (Cc4s::world->rank == 0) {
    elements.push_back(std::make_pair(i, F(1)));
  }
  FockVector<F> basis(
    std::vector<PTR(CTF::Tensor<double>)>(
      {NEW(CTF::Tensor<double>, Tai), NEW(CTF::Tensor<double>, Tabij)}
    ),
    std::vector<std::string>({"ai", "abij"})
  );
  basis *= 0.0;
  basis.write(elements);
  Tai["ai"] = (*basis.get(0))["ai"];
  Tabij["abij"] = (*basis.get(1))["abij"];
}

// instantiate template method implementation
template
void CcsdEquationOfMotionDavidson::getCanonicalPerturbationBasis(
  CTF::Tensor<double> &Tai, CTF::Tensor<double> &Tabij, size_t i
);


template <typename F>
CcsdSimilarityTransformedHamiltonian<F>::CcsdSimilarityTransformedHamiltonian(
  CTF::Tensor<F> *Tai_,
  CTF::Tensor<F> *Tabij_,
  CTF::Tensor<F> *Fij_,
  CTF::Tensor<F> *Fab_,
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
  buildIntermediates();
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

  // Filter out non-physical part
  (*LHijab)["iicd"] = ( 0.0 );
  (*LHijab)["ijcc"] = ( 0.0 );
  (*LHijab)["iicc"] = ( 0.0 );

  return LH;
}

template <typename F>
void CcsdSimilarityTransformedHamiltonian<F>::buildIntermediates() {

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

  auto Tau_abij(NEW(CTF::Tensor<>, *Tabij));
  (*Tau_abij)["abij"] += (*Tai)["ai"] * (*Tai)["bj"];
  (*Tau_abij)["abij"] += ( - 1.0 ) * (*Tai)["bi"] * (*Tai)["aj"];

  //This approach defines intermediates:
  //Wab Wia Wabcd Wabci Waibc
  //Wiabj Wiajk Wij Wijka Wijkl

  Wab   = NEW(CTF::Tensor<>, *Fab);
  Wij   = NEW(CTF::Tensor<>, *Fij);
  Wia   = NEW(CTF::Tensor<>, *Fia);
  Wabcd = NEW(CTF::Tensor<>, *Vabcd);
  Wabci = NEW(CTF::Tensor<>, *Vabci);
  Waibc = NEW(CTF::Tensor<>, *Vaibc);
  Wiabj = NEW(CTF::Tensor<>, *Viabj);
  Wiajk = NEW(CTF::Tensor<>, *Viajk);
  Wijka = NEW(CTF::Tensor<>, *Vijka);
  Wijkl = NEW(CTF::Tensor<>, *Vijkl);

  //Wab
  (*Wab)["ab"]  = (*Fab)["ab"];
  (*Wab)["ab"] += (*Vaibc)["aibc"] * (*Tai)["ci"];
  //(*Wab)["ab"] += ( -1.0) * (*Fia)["ib"] * (*Tai)["ai"];
  (*Wab)["ab"] += (- 0.5) * (*Vijab)["ijbc"] * (*Tabij)["acij"];
  (*Wab)["ab"] += (- 0.5) * (*Vijab)["ijbc"] * (*Tai)["ai"] * (*Tai)["cj"];

  //Wia ( we need this one to construct the 2-body-amplitudes, not directly )
  (*Wia)["ia"] = (*Vijab)["imae"] * (*Tai)["em"];

  //Wij
  (*Wij)["ij"]  = (*Fij)["ij"];
  (*Wij)["ij"] += (*Vijka)["imje"] * (*Tai)["em"];
  //(*Wij)["ij"] += (*Fia)["ie"] * (*Tai)["ej"];
  (*Wij)["ij"] += (  0.5) * (*Vijab)["ikab"] * (*Tabij)["abjk"];
  (*Wij)["ij"] += (  0.5) * (*Vijab)["ikab"] * (*Tai)["ai"] * (*Tai)["bk"];

  //Wabcd
  (*Wabcd)["abcd"]  = (*Vabcd)["abcd"];
  (*Wabcd)["abcd"] += (-1.0) * (*Vaibc)["aicd"] * (*Tai)["bi"];
  // P(ab) for making them symmetric
  (*Wabcd)["abcd"] += ( 1.0) * (*Vaibc)["bicd"] * (*Tai)["ai"];
  (*Wabcd)["abcd"] += ( 0.5) * (*Vijab)["ijcd"] * (*Tai)["ai"] * (*Tai)["bj"];
  (*Wabcd)["abcd"] += ( 0.5) * (*Vijab)["ijcd"] * (*Tabij)["abij"];

  //Wabci TODO: REVIEW
  (*Wabci)["abci"]  = (*Vabci)["abci"];
  (*Wabci)["abci"] += (*Vabcd)["abce"] * (*Tai)["ei"];
  (*Wabci)["abci"] += ( -1.0) * (*Vaicj)["amci"] * (*Tai)["bm"];
  (*Wabci)["abci"] += ( -1.0) * (*Vaicj)["amce"] * (*Tai)["bm"] * (*Tai)["ei"];
  (*Wabci)["abci"] += ( -1.0) * (*Vijak)["mnci"] * (*Tai)["am"] * (*Tai)["bn"];
  //(*Wabci)["abci"] += ( -1.0) * Fia ..... (canonical orbitals)
  (*Wabci)["abci"] += (*Vaibc)["amce"] * (*Tabij)["ebmi"];
  (*Wabci)["abci"] += (*Vijak)["mnci"] * (*Tabij)["abmn"];
  (*Wabci)["abci"] += ( -1.0) * (*Vijab)["mnec"] * (*Tai)["em"] * (*Tabij)["abni"];
  (*Wabci)["abci"] += ( -1.0) * (*Vijab)["mnce"] * (*Tai)["am"] * (*Tabij)["ebni"];
  (*Wabci)["abci"] += (*Vijab)["mnce"] * (*Tai)["ei"] * (*Tabij)["abmn"];
  (*Wabci)["abci"] += (*Vijab)["mnce"] * (*Tai)["am"] * (*Tai)["bn"] * (*Tai)["ei"];

  //Waibc
  (*Waibc)["aibc"]  = (*Vaibc)["aibc"];
  (*Waibc)["aibc"] += ( -1.0) * (*Vijab)["mibc"] * (*Tai)["am"];

  //Wiabj
  (*//)[1] diagram (10.73)
  //This is not listed in the source book, however we can write it in terms
  //of Waijb since it should also have the simmetry of the Tabij amplitudes
  //and the Coulomb integrals Vpqrs
  (*Wiabj)["jabi"]  = (*Vaijb)["ajib"];
  (*Wiabj)["jabi"] += ( -1.0) * (*Vijka)["mjib"] * (*Tai)["am"];
  (*Wiabj)["jabi"] += (*Vaibc)["ajeb"] * (*Tai)["ei"];
  (*Wiabj)["jabi"] += ( -1.0) * (*Vijab)["mjeb"] * (*Tai)["ei"] * (*Tai)["am"];
  (*Wiabj)["jabi"] += (*Vijab)["mjeb"] * (*Tabij)["aeim"];


  //Wijka
  //Taken directly (*from )[2]
  (*Wijka)["ijka"]  = (*Vijka)["jkia"];
  (*Wijka)["ijka"] += (*Tai)["ei"] * (*Vijab)["jkea"];

  //Wijkl
  //Taken directly (*from )[2]
  (*Wijkl)["klij"]  = (*Vijkl)["klij"];
  (*Wijkl)["klij"] += (*Tai)["ej"] * (*Vijka)["klie"];
  (*Wijkl)["klij"] += ( 0.5 ) * (*Tau_abij)["efij"] * (*Vijab)["klef"];


  //Wiajk
  //This is built upon the already existing amplitudes
  (*//)[1] diagram (10.79)
  //Takend directly (*from )[2]
  (*Wiajk)["iajk"]  = (*Viajk)["iajk"];

  //----------------------------------------------------
  (*Wiajk)["iajk"] += (*Vijka)["imje"] * (*Tabij)["aekm"];
  // P(ij)
  (*Wiajk)["iajk"] += ( -1.0 ) * (*Vijka)["jmie"] * (*Tabij)["aekm"];
  //----------------------------------------------------

  (*Wiajk)["iajk"] += (  0.5 ) * (*Viabc)["iaef"] * (*Tau_abij)["efjk"];
  (*Wiajk)["iajk"] += (*Wia)["ie"] * (*Tabij)["aejk"];
  (*Wiajk)["iajk"] += (*Tai)["am"] * (*Vijkl)["imjk"];

  //----------------------------------------------------
  (*Wiajk)["iajk"] += ( -1.0 ) * (*Tai)["ej"] * (*Viabj)["iaek"];
  // P(ij)
  (*Wiajk)["iajk"] += ( +1.0 ) * (*Tai)["ei"] * (*Viabj)["jaek"];
  //----------------------------------------------------

  (*Wiajk)["iajk"] += ( -1.0 ) * (*Tabij)["afmk"] * (*Vijab)["imef"];


}

template <typename F>
FockVector<F> CcsdSimilarityTransformedHamiltonian<F>::rightApply(
  FockVector<F> &R
) {
  FockVector<F> HR(R);
  // get pointers to the component tensors
  PTR(CTF::Tensor<F>) Rai( R.get(0) );
  PTR(CTF::Tensor<F>) Rabij( R.get(1) );
  PTR(CTF::Tensor<F>) HRai( HR.get(0) );
  PTR(CTF::Tensor<F>) HRabij( HR.get(1) );


  // RIGHT APPLY BEGIN
  HRai["ai"]  = 0.0;
  HRai["ai"] += Wab["ad"] * Rai["di"];
  HRai["ai"] += (- 1.0) * Wij["li"] * Rai["di"];
  HRai["ai"] += Wiabj["ladi"] * Rai["dl"];
  HRai["ai"] += Wij["ld"] * Rabij["adil"];
  HRai["ai"] += (   0.5 ) * Waibc["alde"] * Rabij["deil"];
  HRai["ai"] += ( - 0.5 ) * Wijka["lmid"] * Rabij["adlm"];

  // TODO: Deal with permutations automatically
  HRabij["abij"]  = 0.0;

  HRabij["abij"] += Wabci["abej"] * Rai["ei"];
  //P(ij)
  HRabij["abij"] += (-1.0) * Wabci["abei"] * Rai["ej"];

  HRabij["abij"] += (- 1.0 ) * Wiajk["lbij"] * Rai["al"];
  //P(ab)
  HRabij["abij"] += Wiajk["laij"] * Rai["bl"];

  HRabij["abij"] +=           Wab["bd"] * Rabij["adij"];
  //P(ab)
  HRabij["abij"] += ( -1.0) * Wab["ad"] * Rabij["bdij"];

  HRabij["abij"] += ( -1.0) * Wij["lj"] * Rabij["abil"];
  //P(ij)
  HRabij["abij"] += Wij["li"] * Rabij["abjl"];

  HRabij["abij"] += (  0.5) * Wabcd["abde"] * Rabij["deij"];
  HRabij["abij"] += (  0.5) * Wijkl["lmij"] * Rabij["ablm"];

  HRabij["abij"] +=            Wiabj["lbdj"] * Rabij["adil"];
  //-P(ij)
  HRabij["abij"] +=  ( -1.0) * Wiabj["lbdi"] * Rabij["adjl"];
  //-P(ab)
  HRabij["abij"] +=  ( -1.0) * Wiabj["ladj"] * Rabij["bdil"];
  //P(ij)P(ab)
  HRabij["abij"] +=            Wiabj["ladi"] * Rabij["bdjl"];

  // Three body terms from [1]
  // They would have been to be evaluated if we apply it strictly
  //
  //    HRabij["abij"] += Wiabcjk["mabeij"] * Rai["em"];
  //    HRabij["abij"] += (  0.5) *  Waibcdj["ambfej"] * Rabij["feim"];
  //    //P(ij)
  //    HRabij["abij"] += ( -0.5) *  Waibcdj["ambfei"] * Rabij["fejm"];
  //    HRabij["abij"] += ( -0.5) *  Wijakbl["lmbidj"] * Rabij["adlm"];
  //    //P(ab)
  //    HRabij["abij"] += ( 0.5) * Wijakbl["lmaidj"] * Rabij["bdlm"];
  //
  //From [2] equation (31) we can get however another representation
  //involving the T amplitudes, the only part in the right apply
  //that involves the T amplitudes.
  HRabij["abij"] += Rai["em"] * Vaibc["bmfe"] * Tabij["afij"];
  // P(ab)
  HRabij["abij"] += ( -1.0) * Rai["em"] * Vaibc["amfe"] * Tabij["bfij"];

  HRabij["abij"] += ( -0.5) * Rabij["eamn"] * Vijab["nmfe"] * Tabij["fbij"];
  // P(ab)
  HRabij["abij"] += ( -0.5) * Rabij["ebmn"] * Vijab["nmfe"] * Tabij["faij"];

  HRabij["abij"] += ( -1.0) * Rai["em"] * Vijka["nmje"] * Tabij["abin"];
  // P(ij)
  HRabij["abij"] += ( +1.0) * Rai["em"] * Vijka["nmie"] * Tabij["abjn"];

  HRabij["abij"] += ( +0.5) * Rabij["feim"] * Vijab["nmfe"] * Tabij["abjn"];
  // P(ij)
  HRabij["abij"] += ( -0.5) * Rabij["fejm"] * Vijab["nmfe"] * Tabij["abin"];

  // Filter out non-physical part
  //(*HRabij)["cdii"] = ( 0.0 );
  //(*HRabij)["ccij"] = ( 0.0 );
  //(*HRabij)["ccii"] = ( 0.0 );

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
  CTF::Tensor<F> &Vijkl
): diagonalH(
    std::vector<PTR(CTF::Tensor<double>)>(
      {NEW(CTF::Tensor<double>, Tai), NEW(CTF::Tensor<double>, Tabij)}
    ),
    std::vector<std::string>({"ai", "abij"})
  ) {
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

  // Filter out non-physical part
  (*Dabij)["cdii"] = ( 0.0 );
  (*Dabij)["ccij"] = ( 0.0 );
  (*Dabij)["ccii"] = ( 0.0 );
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
std::vector<FockVector<F>> CcsdPreConditioner<F>::getInitialBasis(
  const int eigenVectorsCount
) {
  LOG(0, "CcsdEomDavid") << "Getting initial basis " << std::endl;
  // find K=eigenVectorsCount lowest diagonal elements at each processor
  std::vector<std::pair<size_t, F>> localElements( diagonalH.readLocal() );
  std::sort(
    localElements.begin(), localElements.end(),
    EomDiagonalValueComparator<double>()
  );

  // gather all K elements of all processors at root
  //   convert into homogeneous arrays for MPI gather
  std::vector<size_t> localLowestElementIndices(localElements.size());
  std::vector<F> localLowestElementValues(localElements.size());
  for (size_t i(0); i < localElements.size(); ++i) {
    localLowestElementIndices[i] = localElements[i].first;
    localLowestElementValues[i] = localElements[i].second;
  }
  MpiCommunicator communicator(*Cc4s::world);
   int lowestElementsCount(
    diagonalH.get(0)->lens[0] *
    diagonalH.get(0)->lens[1] +
    pow(
      diagonalH.get(0)->lens[0] *
      diagonalH.get(0)->lens[1],
      3.0
    )
  );
  std::vector<size_t> lowestElementIndices(lowestElementsCount);
  std::vector<F> lowestElementValues(lowestElementsCount);
  communicator.gather(localLowestElementIndices, lowestElementIndices);
  communicator.gather(localLowestElementValues, lowestElementValues);
  //   convert back into (index,value) pairs for sorting
  std::vector<std::pair<size_t, F>> lowestElements(lowestElementsCount);
  for (int i(0); i < lowestElementsCount; ++i) {
    lowestElements[i].first = lowestElementIndices[i];
    lowestElements[i].second = lowestElementValues[i];
  }

  // find globally lowest K diagonal elements among the gathered elements
  std::sort(
    lowestElements.begin(), lowestElements.end(),
    EomDiagonalValueComparator<double>()
  );
  // at rank==0 (root) lowestElements contains N*Np entries
  // rank > 0 has an empty list

  // create basis vectors for each lowest element
  std::vector<V> basis;
  //for (int b(0); b < eigenVectorsCount; ++b) {
  int bb(0);
  int b(0);
  while (bb < eigenVectorsCount) {
    V basisElement(diagonalH);
    basisElement *= 0.0;
    std::vector<std::pair<size_t,F>> elements;
    if (communicator.getRank() == 0) {
      elements.push_back(
        std::make_pair(lowestElements[b].first, 1.0)
      );
    }
    basisElement.write(elements);
    // (101, -70), (32, -55), ...
    // b1: 0... 1 (at global position 101) 0 ...
    // b2: 0... 1 (at global position 32) 0 ...i

    // Filter out unphysical components from the basisElement
    (*basisElement.get(1))["abii"]=0.0;
    (*basisElement.get(1))["aaij"]=0.0;
    (*basisElement.get(1))["aaii"]=0.0;

    b++;
    //std::cout << "b" << b << std::endl;
    if (std::sqrt(basisElement.dot(basisElement))!=F(1)) continue;
    bb++;
    basis.push_back(basisElement);
    //std::cout << "bb" << bb << std::endl;
  }
  return basis;
}

template <typename F>
FockVector<F> CcsdPreConditioner<F>::getCorrection(
  const complex lambda, FockVector<F> &residuum
) {
  FockVector<F> w(diagonalH);

  // Define a helping class for the diagonal correction
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
  return correction;
}


// instantiate class
template
class CcsdPreConditioner<double>;

