#include <algorithms/Mp2EquationOfMotionDavidson.hpp>

#include <tcc/Tcc.hpp>
#include <tcc/DryMachineTensor.hpp>
#include <math/EigenSystemDavidson.hpp>
#include <math/MathFunctions.hpp>
#include <math/FockVector.hpp>
#include <math/RandomTensor.hpp>
#include <util/MpiCommunicator.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>

#include <algorithm>
#include <utility>

using namespace cc4s;
using namespace tcc;

ALGORITHM_REGISTRAR_DEFINITION(Mp2EquationOfMotionDavidson);

Mp2EquationOfMotionDavidson::Mp2EquationOfMotionDavidson(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}
Mp2EquationOfMotionDavidson::~Mp2EquationOfMotionDavidson() {}

void Mp2EquationOfMotionDavidson::run() {
  typedef CTF::Tensor<> T;

  // Get copy of couloumb integrals
  T Vijkl(*getTensorArgument<double, T>("HHHHCoulombIntegrals"));
  T Vabcd(*getTensorArgument<double, T>("PPPPCoulombIntegrals")); 

  T Vabij(*getTensorArgument<double, T>("PPHHCoulombIntegrals"));
  // T *Vijab(getTensorArgument<double, T>("HHPPCoulombIntegrals")); // swap PPHH (done)

  T Vijka(*getTensorArgument<double, T>("HHHPCoulombIntegrals"));
  // T *Viajk(getTensorArgument<double, T>("HPHHCoulombIntegrals")); // swap HHHP (done)

  T Vaibj(*getTensorArgument<double, T>("PHPHCoulombIntegrals")); // not in eqs
  //T *Viajb(getTensorArgument<double, T>("HPHPCoulombIntegrals")); // swap PHPH (done)

  T Vabci(*getTensorArgument<double, T>("PPPHCoulombIntegrals")); // not in eqs
  //T *Viabc(getTensorArgument<double, T>("HPPPCoulombIntegrals")); // swap PPPH (done)
  //T *Vabic(getTensorArgument<double, T>("PPHPCoulombIntegrals")); // swap PPPH (done)


  // Get orbital energies
  T *epsi(getTensorArgument<double, T>("HoleEigenEnergies"));
  T *epsa(getTensorArgument<double, T>("ParticleEigenEnergies"));
  int Nv(epsa->lens[0]), No(epsi->lens[0]);
  int syms[] = {NS, NS, NS, NS};

  //  Tai
  int vo[] = { Nv, No };
  T Tai(2, vo, syms, *Cc4s::world, "Tai");
  T Tabij(false, Vabij);

  LOG(1, "MP2_EOM_DAVIDSON") << "Antisymmetrizing Vpqrs" << std::endl;

  //  Vijab
  int oovv[] = { No, No, Nv, Nv };
  T Vijab(4, oovv, syms, *Cc4s::world, "Vijab");
  Vijab["ijab"] =  Vabij["abij"] - Vabij["abji"];

  //  Viajk
  int ovoo[] = { No, Nv, No, No };
  T Viajk(4, ovoo, syms, *Cc4s::world, "Viajk");
  Viajk["iajk"] =  Vijka["ijka"]  - Vijka["ikja"];

  // Viajb
  int ovov[] = { No, Nv, No, Nv };
  T Viajb(4, ovov, syms, *Cc4s::world, "Viajb");
  Viajb["iajb"] =  Vaibj["aibj"] - Vabij["abji"];

  // Viabc
  int ovvv[] = { No, Nv, Nv, Nv };
  T Viabc(4, ovvv, syms, *Cc4s::world, "Viabc");
  Viabc["iabc"] =  Vabci["abci"] - Vabci["acbi"];

  // Vabic
  int vvov[] = { Nv, Nv, No, Nv };
  T Vabic(4, vvov, syms, *Cc4s::world, "Vabic");
  Vabic["abic"] =  Vabci["abci"] - Vabci["baci"];

  // Antisymmetrize integrals
  Vijkl["ijkl"] -= Vijkl["ijlk"];
  Vabcd["abcd"] -= Vabcd["abdc"];
  Vijka["ijka"] -= Vijka["jika"];
  Vaibj["aibj"] -= Vabij["baij"];
  Vabci["abci"] -= Vabci["baci"];
  Vabij["abij"] -= Vabij["abji"];

  Tabij["abij"] =  (*epsi)["i"];
  Tabij["abij"] += (*epsi)["j"];
  Tabij["abij"] -= (*epsa)["a"];
  Tabij["abij"] -= (*epsa)["b"];

  LOG(1, "MP2_EOM_DAVIDSON") << "Creating doubles amplitudes" << std::endl;
  CTF::Bivar_Function<> fDivide(&divide<double>);
  Tabij.contract(1.0, Vabij,"abij", Tabij,"abij", 0.0,"abij", fDivide);

  CTF::Scalar<> energy(0.0);
  double e(0.0);
  LOG(2, "MP2_EOM_DAVIDSON") << "Calculating MP2 energy" << std::endl;
  energy[""] = ( 0.25 ) * Tabij["abij"] * Vabij["abij"];
  e = energy.get_val();
  LOG(1, "MP2_EOM_DAVIDSON") << " Mp2 energy = " << e << std::endl;

  // HF terms
  int kineticLensVirtual[] = {Nv, Nv};
  int kineticSyms[] = {NS, NS};
  T Fab(2, kineticLensVirtual, kineticSyms, *Cc4s::world, "Fab");
  int kineticLensOccupied[] = {No, No};
  T Fij(2, kineticLensOccupied, kineticSyms, *Cc4s::world, "Fij");
  Fab["aa"] = (*epsa)["a"];
  Fij["ii"] = (*epsi)["i"];

  Mp2SimilarityTransformedHamiltonian<double> H(
    &Tai, &Tabij, &Fij, &Fab,
    &Vabcd, &Viajb, &Vijab, &Vijkl, &Vijka, &Viabc, &Viajk, &Vabic
    // TODO add extra Vs...
  );
  Mp2PreConditioner<double> P(
    Tai, Tabij, Fij, Fab, Vabcd, Viajb, Vijab, Vijkl
  );

  // .. davidson solver, needing Mp2PreConditioner & Mp2SimTransHam

}

// template method implementation
template <typename F>
void Mp2EquationOfMotionDavidson::getCanonicalPerturbationBasis(
  CTF::Tensor<F> &Tai, CTF::Tensor<F> &Tabij, int64_t i
) {
  std::vector<std::pair<int64_t, F>> elements;
  if (Cc4s::world->rank == 0) {
    elements.push_back(std::make_pair(i, F(1)));
  }
  FockVector<F> basis({Tai, Tabij}, {"ai", "abij"});
  basis *= 0.0;
  basis.write(elements);
  Tai["ai"] = basis.componentTensors[0]["ai"];
  Tabij["abij"] = basis.componentTensors[1]["abij"];
}

// instantiate template method implementation
template
void Mp2EquationOfMotionDavidson::getCanonicalPerturbationBasis(
  CTF::Tensor<double> &Tai, CTF::Tensor<double> &Tabij, int64_t i
);


template <typename F>
Mp2SimilarityTransformedHamiltonian<F>::Mp2SimilarityTransformedHamiltonian(
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
FockVector<F> Mp2SimilarityTransformedHamiltonian<F>::rightApply(
  FockVector<F> &R
) {
  FockVector<F> HR(R);
  // get pointers to the component tensors
  CTF::Tensor<F> *Rai( &R.componentTensors[0] );
  CTF::Tensor<F> *Rabij( &R.componentTensors[1] );
  CTF::Tensor<F> *HRai( &HR.componentTensors[0] );
  CTF::Tensor<F> *HRabij( &HR.componentTensors[1] );

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

  // symmetrize doubles part w.r.t. left/right
  (*HRabij)["abij"] += (*HRabij)["baji"];
  // anti-symmetrize w.r.t. exchange
  (*HRabij)["abij"] -= (*HRabij)["abji"];
  (*HRabij)["abij"] *= 0.5;

  return HR;
}

// instantiate class
template
class Mp2SimilarityTransformedHamiltonian<double>;


template <typename F>
Mp2PreConditioner<F>::Mp2PreConditioner(
  CTF::Tensor<F> &Tai,
  CTF::Tensor<F> &Tabij,
  CTF::Tensor<F> &Fij,
  CTF::Tensor<F> &Fab,
  CTF::Tensor<F> &Vabcd,
  CTF::Tensor<F> &Viajb,
  CTF::Tensor<F> &Vijab,
  CTF::Tensor<F> &Vijkl
): diagonalH({{Tai, Tabij}, {"ai", "abij"}}) {
  // pointers to singles and doubles tensors of diagonal part
  CTF::Tensor<F> *Dai( &diagonalH.componentTensors[0] );
  CTF::Tensor<F> *Dabij( &diagonalH.componentTensors[1] );

  // calculate diagonal elements of H
  (*Dai)["aj"]  = ( - 1.0  ) * Fij["jj"];
  (*Dai)["aj"] += ( + 1.0  ) * Fab["aa"];
  (*Dai)["aj"] += ( - 1.0  ) * Viajb["jaja"];
  (*Dai)["aj"] += ( + 1.0  ) * Tabij["adjn"] * Vijab["njda"];
  (*Dai)["aj"] += ( + 0.5  ) * Tabij["cdjn"] * Vijab["njcd"];
  (*Dai)["aj"] += ( + 0.5  ) * Tabij["admn"] * Vijab["mnda"];

  (*Dabij)["abkl"] += ( - 1.0  ) * Fij["kk"];
/*
  (*Dabij)["abkl"] += ( + 1.0  ) * Rabij["abkl"] * Fij["lm"] * Lijab["mkab"];
  (*Dabij)["abkl"] += ( - 1.0  ) * Rabij["abkl"] * Fab["eb"] * Lijab["klea"];
  (*Dabij)["abkl"] += ( + 1.0  ) * Rabij["abkl"] * Fab["ea"] * Lijab["kleb"];
  (*Dabij)["abkl"] += ( - 0.5  ) * Rabij["abkl"] * Vijkl["klmn"] * Lijab["nmab"];
  (*Dabij)["abkl"] += ( + 1.0  ) * Rabij["abkl"] * Viajb["kenb"] * Lijab["nlea"];
  (*Dabij)["abkl"] += ( - 1.0  ) * Rabij["abkl"] * Viajb["kena"] * Lijab["nleb"];
  (*Dabij)["abkl"] += ( - 1.0  ) * Rabij["abkl"] * Viajb["lenb"] * Lijab["nkea"];
  (*Dabij)["abkl"] += ( + 1.0  ) * Rabij["abkl"] * Viajb["lena"] * Lijab["nkeb"];
  (*Dabij)["abkl"] += ( - 0.5  ) * Rabij["abkl"] * Vabcd["efab"] * Lijab["klfe"];
  (*Dabij)["abkl"] += ( + 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["klfb"] * Lijab["poea"];
  (*Dabij)["abkl"] += ( - 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["klfa"] * Lijab["poeb"];
  (*Dabij)["abkl"] += ( - 0.25  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["klef"] * Lijab["poab"];
  (*Dabij)["abkl"] += ( - 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["pkab"] * Lijab["olfe"];
  (*Dabij)["abkl"] += ( + 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["plab"] * Lijab["okfe"];
  (*Dabij)["abkl"] += ( - 1.0  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["pkfb"] * Lijab["olea"];
  (*Dabij)["abkl"] += ( + 1.0  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["pkfa"] * Lijab["oleb"];
  (*Dabij)["abkl"] += ( + 1.0  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["plfb"] * Lijab["okea"];
  (*Dabij)["abkl"] += ( - 1.0  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["plfa"] * Lijab["okeb"];
  (*Dabij)["abkl"] += ( + 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["pkef"] * Lijab["olab"];
  (*Dabij)["abkl"] += ( - 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["plef"] * Lijab["okab"];
  (*Dabij)["abkl"] += ( - 0.25  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["opab"] * Lijab["klfe"];
  (*Dabij)["abkl"] += ( - 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["opfb"] * Lijab["klea"];
  (*Dabij)["abkl"] += ( + 0.5  ) * Rabij["abkl"] * Tabij["efop"] * Vijab["opfa"] * Lijab["kleb"];
*/
  // TODO: antisymmetrize
}


template <typename F>
std::vector<FockVector<F>> Mp2PreConditioner<F>::getInitialBasis(
  const int eigenVectorsCount
) {
  // find K=eigenVectorsCount lowest diagonal elements at each processor
  std::vector<std::pair<int64_t, F>> localElements( diagonalH.readLocal() );
  std::sort(
    localElements.begin(), localElements.end(),
    typename LapackGeneralEigenSystem<F>::EigenValueComparator()
  );

  // gather all K elements of all processors at root
  //   convert into homogeneous arrays for MPI gather
  std::vector<int64_t> localLowestElementIndices;
  std::vector<F> localLowestElementValues;
  for (int64_t i(0); i < eigenVectorsCount; ++i) {
    localLowestElementIndices[i] = localElements[i].first;
    localLowestElementValues[i] = localElements[i].second;
  }
  MpiCommunicator communicator(*Cc4s::world);
  int lowestElementsCount(
    communicator.getRank() == 0 ?
      eigenVectorsCount * communicator.getProcesses() : 0
  );
  std::vector<int64_t> lowestElementIndices(lowestElementsCount);
  std::vector<F> lowestElementValues(lowestElementsCount);
  communicator.gather(localLowestElementIndices, lowestElementIndices);
  communicator.gather(localLowestElementValues, lowestElementValues);
  //   convert back into (index,value) pairs for sorting
  std::vector<std::pair<int64_t, F>> lowestElements(lowestElementsCount);
  for (int i(0); i < lowestElementsCount; ++i) {
    lowestElements[i].first = lowestElementIndices[i];
    lowestElements[i].second = lowestElementValues[i];
  }

  // find globally lowest K diagonal elements among the gathered elements
  std::sort(
    lowestElements.begin(), lowestElements.end(),
    typename LapackGeneralEigenSystem<F>::EigenValueComparator()
  );
  // at rank==0 (root) lowestElements contains N*Np entries
  // rank > 0 has an empty list

  // create basis vectors for each lowest element
  std::vector<V> basis;
  for (int b(0); b < eigenVectorsCount; ++b) {
    V basisElement(diagonalH);
    basisElement *= 0.0;
    std::vector<std::pair<int64_t,F>> elements;
    if (communicator.getRank() == 0) {
      elements.push_back(
        std::make_pair(lowestElements[b].first, 1.0)
      );
    }
    basisElement.write(elements);
    // (101, -70), (32, -55), ...
    // b1: 0... 1 (at global position 101) 0 ...
    // b2: 0... 1 (at global position 32) 0 ...i
    basis.push_back(basisElement);
  }
  return basis;
}

template <typename F>
FockVector<F> Mp2PreConditioner<F>::getCorrection(
  const complex lambda, FockVector<F> &residuum
) {
  FockVector<F> w(diagonalH);
  class DiagonalCorrection {
  public:
    DiagonalCorrection(const double lambda_): lambda(lambda_) {
    }
    F operator ()(const F residuumElement, const F diagonalElement) {
      return residuumElement / (lambda - diagonalElement);
    }
  protected:
    double lambda;
  } diagonalCorrection(std::real(lambda));

  FockVector<F> correction(diagonalH);
  // compute ((lambda * id - Diag(diagonal))^-1) . residuum
  for (unsigned int c(0); c < w.componentTensors.size(); ++c) {
    const char *indices( correction.componentIndices[c].c_str() );
    correction.componentTensors[c].contract(
      1.0,
      residuum.componentTensors[c],indices,
      diagonalH.componentTensors[c],indices,
      0.0,indices,
      CTF::Bivar_Function<F>(diagonalCorrection)
    );
  }
  return correction;
}


// instantiate class
template
class Mp2PreConditioner<double>;

