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

  //CTF::Tensor<> *Vaibj(
      //getTensorArgument<double, CTF::Tensor<>>("PHPHCoulombIntegrals"));
  //CTF::Tensor<> *Vabij(
      //getTensorArgument<double, CTF::Tensor<>>("PPHHCoulombIntegrals"));
  //CTF::Tensor<> *Vabci(
      //getTensorArgument<double, CTF::Tensor<>>("PPPHCoulombIntegrals"));

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
    Vabcd, Viajb, Vijab, Vijkl, Vijka, Viabc, Viajk, Vabic
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
  LOG(0, "CCSD_EOM_DAVIDSON") << "Computing " << eigenStates << " eigen states"
                              << std::endl;
  EigenSystemDavidson<FockVector<double>> eigenSystem(H, eigenStates, P, 1E-4, 8*16);

  std::vector<complex> eigenValues(eigenSystem.getEigenValues());
  for (auto &ev: eigenValues) {
    LOG(0, "CCSD_EOM_DAVIDSON") << "Eigenvalue=" << ev << std::endl;
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
  CTF::Tensor<F> *Vabic_
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
  Vabic(Vabic_)
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

  // Filter out non-physical part
  (*LHijab)["iicd"] = ( 0.0 );
  (*LHijab)["ijcc"] = ( 0.0 );
  (*LHijab)["iicc"] = ( 0.0 );

  return LH;
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

  // Contruct HR (one body part)
  // TODO: why "bi" not "ai"?
  (*HRai)["bi"]  = 0.0;
  (*HRai)["bi"] += ( - 1.0 ) * (*Fij)["ki"] * (*Rai)["bk"];
  (*HRai)["bi"] += ( + 1.0 ) * (*Fab)["bc"] * (*Rai)["ci"];
  (*HRai)["bi"] += ( - 1.0 ) * (*Viajb)["kbid"] * (*Rai)["dk"];
  (*HRai)["bi"] += ( + 0.5 ) * (*Vijka)["klie"] * (*Rabij)["ebkl"];
  (*HRai)["bi"] += ( + 0.5 ) * (*Viabc)["kbde"] * (*Rabij)["deki"];
  (*HRai)["bi"] += ( + 1.0 ) * (*Tabij)["cbli"] * (*Vijab)["lmcf"] * (*Rai)["fm"];
  (*HRai)["bi"] += ( - 0.5 ) * (*Tabij)["cdmi"] * (*Vijab)["mncd"] * (*Rai)["bn"];
  (*HRai)["bi"] += ( - 0.5 ) * (*Tabij)["cblm"] * (*Vijab)["lmcf"] * (*Rai)["fi"];

  // Singles amplitudes part
  (*HRai)["bi"] += ( - 1.0  ) * (*Tai)["bk"] * (*Vijka)["klie"] * (*Rai)["el"];
  (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["cl"] * (*Vijka)["lmic"] * (*Rai)["bm"];
  (*HRai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Viabc)["lbce"] * (*Rai)["el"];
  (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["cl"] * (*Viabc)["lbce"] * (*Rai)["ei"];
  (*HRai)["bi"] += ( + 0.5  ) * (*Tai)["ci"] * (*Vijab)["lmcf"] * (*Rabij)["fblm"];
  (*HRai)["bi"] += ( + 0.5  ) * (*Tai)["bk"] * (*Vijab)["klef"] * (*Rabij)["efli"];
  (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["cl"] * (*Vijab)["lmcf"] * (*Rabij)["fbmi"];
  (*HRai)["bi"] += ( - 1.0  ) * (*Tai)["ci"] * (*Tai)["bl"] * (*Vijab)["lmcf"] * (*Rai)["fm"];
  (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["ci"] * (*Tai)["dm"] * (*Vijab)["mncd"] * (*Rai)["bn"];
  (*HRai)["bi"] += ( + 1.0  ) * (*Tai)["bk"] * (*Tai)["dm"] * (*Vijab)["kmdf"] * (*Rai)["fi"];

  // Contruct HR (two body part)
  // TODO: why "cdai" not "abij"?
  (*HRabij)["cdij"]  = 0.0;
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Viajk)["mdij"] * (*Rai)["cm"];
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Viajk)["mcij"] * (*Rai)["dm"];
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Vabic)["cdie"] * (*Rai)["ej"];
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Vabic)["cdje"] * (*Rai)["ei"];
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Fij)["mi"] * (*Rabij)["cdmj"];
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Fij)["mj"] * (*Rabij)["cdmi"];
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Fab)["de"] * (*Rabij)["ecij"];
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Fab)["ce"] * (*Rabij)["edij"];
  (*HRabij)["cdij"] += ( + 0.5 ) * (*Vijkl)["mnij"] * (*Rabij)["cdmn"];
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Viajb)["mdif"] * (*Rabij)["fcmj"];
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Viajb)["mcif"] * (*Rabij)["fdmj"];
  (*HRabij)["cdij"] += ( - 1.0 ) * (*Viajb)["mdjf"] * (*Rabij)["fcmi"];
  (*HRabij)["cdij"] += ( + 1.0 ) * (*Viajb)["mcjf"] * (*Rabij)["fdmi"];
  (*HRabij)["cdij"] += ( + 0.5 ) * (*Vabcd)["cdef"] * (*Rabij)["efij"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["cdmj"] * (*Vijka)["mnig"] * (*Rai)["gn"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["cdmi"] * (*Vijka)["mnjg"] * (*Rai)["gn"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ednj"] * (*Vijka)["noie"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["ecnj"] * (*Vijka)["noie"] * (*Rai)["do"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["edni"] * (*Vijka)["noje"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ecni"] * (*Vijka)["noje"] * (*Rai)["do"];
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["cdmn"] * (*Vijka)["mnig"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["cdmn"] * (*Vijka)["mnjg"] * (*Rai)["gi"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ecij"] * (*Viabc)["ndeg"] * (*Rai)["gn"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["edij"] * (*Viabc)["nceg"] * (*Rai)["gn"];
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["efij"] * (*Viabc)["odef"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["efij"] * (*Viabc)["ocef"] * (*Rai)["do"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ecni"] * (*Viabc)["ndeg"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["edni"] * (*Viabc)["nceg"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["ecnj"] * (*Viabc)["ndeg"] * (*Rai)["gi"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ednj"] * (*Viabc)["nceg"] * (*Rai)["gi"];
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["edij"] * (*Vijab)["noeh"] * (*Rabij)["hcno"];
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["ecij"] * (*Vijab)["noeh"] * (*Rabij)["hdno"];
  (*HRabij)["cdij"] +=
    ( + 0.25) * (*Tabij)["efij"] * (*Vijab)["opef"] * (*Rabij)["cdop"];
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["cdmi"] * (*Vijab)["mngh"] * (*Rabij)["ghnj"];
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["cdmj"] * (*Vijab)["mngh"] * (*Rabij)["ghni"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["edni"] * (*Vijab)["noeh"] * (*Rabij)["hcoj"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ecni"] * (*Vijab)["noeh"] * (*Rabij)["hdoj"];
  (*HRabij)["cdij"] +=
    ( + 1.0 ) * (*Tabij)["ednj"] * (*Vijab)["noeh"] * (*Rabij)["hcoi"];
  (*HRabij)["cdij"] +=
    ( - 1.0 ) * (*Tabij)["ecnj"] * (*Vijab)["noeh"] * (*Rabij)["hdoi"];
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["efoi"] * (*Vijab)["opef"] * (*Rabij)["cdpj"];
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["efoj"] * (*Vijab)["opef"] * (*Rabij)["cdpi"];
  (*HRabij)["cdij"] +=
    ( + 0.25) * (*Tabij)["cdmn"] * (*Vijab)["mngh"] * (*Rabij)["ghij"];
  (*HRabij)["cdij"] +=
    ( + 0.5 ) * (*Tabij)["edno"] * (*Vijab)["noeh"] * (*Rabij)["hcij"];
  (*HRabij)["cdij"] +=
    ( - 0.5 ) * (*Tabij)["ecno"] * (*Vijab)["noeh"] * (*Rabij)["hdij"];

  // Singles amplitudes part
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["dm"] * (*Vijkl)["mnij"] * (*Rai)["cn"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["cm"] * (*Vijkl)["mnij"] * (*Rai)["dn"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Viajb)["ndie"] * (*Rai)["cn"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ej"] * (*Viajb)["ncie"] * (*Rai)["dn"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Viajb)["ndje"] * (*Rai)["cn"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Viajb)["ncje"] * (*Rai)["dn"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["cm"] * (*Viajb)["mdif"] * (*Rai)["fj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["dm"] * (*Viajb)["mcif"] * (*Rai)["fj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["cm"] * (*Viajb)["mdjf"] * (*Rai)["fi"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["dm"] * (*Viajb)["mcjf"] * (*Rai)["fi"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Vabcd)["cdef"] * (*Rai)["fj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Vabcd)["cdef"] * (*Rai)["fi"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tai)["ej"] * (*Vijka)["noie"] * (*Rabij)["cdno"];
  (*HRabij)["cdij"] +=
    ( - 0.5  ) * (*Tai)["ei"] * (*Vijka)["noje"] * (*Rabij)["cdno"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["dm"] * (*Vijka)["mnig"] * (*Rabij)["gcnj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["cm"] * (*Vijka)["mnig"] * (*Rabij)["gdnj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["dm"] * (*Vijka)["mnjg"] * (*Rabij)["gcni"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["cm"] * (*Vijka)["mnjg"] * (*Rabij)["gdni"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["en"] * (*Vijka)["noie"] * (*Rabij)["cdoj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["en"] * (*Vijka)["noje"] * (*Rabij)["cdoi"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Viabc)["ndeg"] * (*Rabij)["gcnj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Viabc)["nceg"] * (*Rabij)["gdnj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Viabc)["ndeg"] * (*Rabij)["gcni"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ej"] * (*Viabc)["nceg"] * (*Rabij)["gdni"];
  (*HRabij)["cdij"] +=
    ( - 0.5  ) * (*Tai)["cm"] * (*Viabc)["mdfg"] * (*Rabij)["fgij"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tai)["dm"] * (*Viabc)["mcfg"] * (*Rabij)["fgij"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["en"] * (*Viabc)["ndeg"] * (*Rabij)["gcij"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["en"] * (*Viabc)["nceg"] * (*Rabij)["gdij"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Vijka)["noie"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Vijka)["noie"] * (*Rai)["do"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Vijka)["noje"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Vijka)["noje"] * (*Rai)["do"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijka)["mnig"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijka)["mnjg"] * (*Rai)["gi"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Viabc)["odef"] * (*Rai)["co"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Viabc)["ocef"] * (*Rai)["do"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Viabc)["ndeg"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Viabc)["nceg"] * (*Rai)["gj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Viabc)["ndeg"] * (*Rai)["gi"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Viabc)["nceg"] * (*Rai)["gi"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Vijab)["opef"] * (*Rabij)["cdop"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["dn"] * (*Vijab)["noeh"] * (*Rabij)["hcoj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Vijab)["noeh"] * (*Rabij)["hdoj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["dn"] * (*Vijab)["noeh"] * (*Rabij)["hcoi"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Vijab)["noeh"] * (*Rabij)["hdoi"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fo"] * (*Vijab)["opef"] * (*Rabij)["cdpj"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["fo"] * (*Vijab)["opef"] * (*Rabij)["cdpi"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tai)["cm"] * (*Tai)["dn"] * (*Vijab)["mngh"] * (*Rabij)["ghij"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["dm"] * (*Tai)["fo"] * (*Vijab)["mofh"] * (*Rabij)["hcij"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["cm"] * (*Tai)["fo"] * (*Vijab)["mofh"] * (*Rabij)["hdij"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["cdmj"] * (*Tai)["fi"] * (*Vijab)["mofh"] * (*Rai)["ho"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["cdmi"] * (*Tai)["fj"] * (*Vijab)["mofh"] * (*Rai)["ho"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["ednj"] * (*Tai)["gi"] * (*Vijab)["npeg"] * (*Rai)["cp"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["gi"] * (*Vijab)["npeg"] * (*Rai)["dp"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["edni"] * (*Tai)["gj"] * (*Vijab)["npeg"] * (*Rai)["cp"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["ecni"] * (*Tai)["gj"] * (*Vijab)["npeg"] * (*Rai)["dp"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gi"] * (*Vijab)["mngh"] * (*Rai)["hj"];
  (*HRabij)["cdij"] +=
    ( - 0.5  ) * (*Tabij)["cdmn"] * (*Tai)["gj"] * (*Vijab)["mngh"] * (*Rai)["hi"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["ecij"] * (*Tai)["dn"] * (*Vijab)["noeh"] * (*Rai)["ho"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["edij"] * (*Tai)["cn"] * (*Vijab)["noeh"] * (*Rai)["ho"];
  (*HRabij)["cdij"] +=
    ( - 0.5  ) * (*Tabij)["efij"] * (*Tai)["do"] * (*Vijab)["opef"] * (*Rai)["cp"];
  (*HRabij)["cdij"] +=
    ( + 0.5  ) * (*Tabij)["efij"] * (*Tai)["co"] * (*Vijab)["opef"] * (*Rai)["dp"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["edij"] * (*Tai)["fo"] * (*Vijab)["opef"] * (*Rai)["cp"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["ecij"] * (*Tai)["fo"] * (*Vijab)["opef"] * (*Rai)["dp"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["ecni"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["edni"] * (*Tai)["co"] * (*Vijab)["noeh"] * (*Rai)["hj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["ecnj"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hi"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["ednj"] * (*Tai)["co"] * (*Vijab)["noeh"] * (*Rai)["hi"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tabij)["cdmi"] * (*Tai)["fo"] * (*Vijab)["mofh"] * (*Rai)["hj"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tabij)["cdmj"] * (*Tai)["fo"] * (*Vijab)["mofh"] * (*Rai)["hi"];
  (*HRabij)["cdij"] +=
    ( - 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["do"] * (*Vijab)["opef"] * (*Rai)["cp"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["fj"] * (*Tai)["co"] * (*Vijab)["opef"] * (*Rai)["dp"];
  (*HRabij)["cdij"] +=
    ( + 1.0  ) * (*Tai)["ei"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hj"];
  (*HRabij)["cdij"] +=
   ( - 1.0  ) * (*Tai)["ej"] * (*Tai)["cn"] * (*Tai)["do"] * (*Vijab)["noeh"] * (*Rai)["hi"];


  // Filter out non-physical part
  (*HRabij)["cdii"] = ( 0.0 );
  (*HRabij)["ccij"] = ( 0.0 );
  (*HRabij)["ccii"] = ( 0.0 );

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
  LOG(0, "CCSD_EOM_DAVIDSON") << "Getting initial basis " << std::endl;
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

