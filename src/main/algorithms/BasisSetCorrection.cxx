#include <algorithms/BasisSetCorrection.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <math/TensorUnion.hpp>


using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(BasisSetCorrection)

Ptr<MapNode> BasisSetCorrection::run(const Ptr<MapNode> &arguments) {
  auto result(New<MapNode>(SOURCE_LOCATION));
  // multiplex calls to template methods
  bool success(false);
  if (Cc4s::options->dryRun) {
    using TE = DefaultDryTensorEngine;
    success =
      run<Real<>,TE>(arguments, result) ||
      run<Complex<>,TE>(arguments, result);
  } else {
    using TE = DefaultTensorEngine;
    success =
      run<Real<>,TE>(arguments, result) ||
      run<Complex<>,TE>(arguments, result);
  }
  ASSERT(
    success, "unsupported orbitals type in amplitudes"
  );
  return result;
}



template <typename F, typename TE>
bool BasisSetCorrection::run(
  const Ptr<MapNode> &arguments, Ptr<MapNode> &result
) {
  using TRc = TensorRecipe<Complex<>, TE>;
  using Tc = Tensor<Complex<>, TE>;
  using Tr = Tensor<Real<>, TE>;
  using T  = Tensor<F, TE>;
  auto amplitudesNode(
    arguments->get("amplitudes")->toAtom<Ptr<const TensorUnion<F,TE>>>()
  );
  //amplitudesNode is nullptr if the node was of different type
  if (!amplitudesNode) return false;
  auto amplitudes(amplitudesNode->value);
  auto Tph( amplitudes->get(0) );
  auto Tpphh( amplitudes->get(1) );

  auto nij(arguments->getValue<Ptr<T>>("nij"));
  auto Dabij(arguments->getValue<Ptr<T>>("deltaIntegrals"));
  // these should be the cbs estimates for the mp2 pair energies
  auto mp2PairEnergiesNode(arguments->getMap("mp2PairEnergies"));
  auto mp2PairEnergiesCbs(mp2PairEnergiesNode->getValue<Ptr<Tr>>("data"));

  // these are the mp2 pair energies in the fno/finite/cc4s basis
  auto mp2PairEnergiesFno( Tcc<TE>::template tensor<F>("mp2PairEnergiesFno"));
  auto eMp2( Tcc<TE>::template tensor<F>("Emp2"));
  auto eMp2Cbs( Tcc<TE>:: template tensor<Real<>>("Emp2Cbs"));

  auto eCcsd( Tcc<TE>::template tensor<F>("Eccsd"));


  auto Tabij( Tcc<TE>::template tensor<F>("Tabij"));

  auto gijccd( Tcc<TE>::template tensor<F>("gijccd"));
  auto gijmp2( Tcc<TE>::template tensor<F>("gijmp2"));
  auto geff  ( Tcc<TE>::template tensor<Real<>>("geff"));
  auto deltaEppl( Tcc<TE>::template tensor<Real<>>("deltaEppl"));
  //mp2 amplitudes
  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto coulombSlices(coulombIntegrals->getMap("slices"));
  auto Vabij(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("pphh"));

  auto eigenEnergies(arguments->getMap("slicedEigenEnergies"));
  auto energySlices(eigenEnergies->getMap("slices"));
  auto epsi(energySlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("h"));
  auto epsa(energySlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("p"));

  auto No(epsi->getResult()->lens[0]);
  auto Nv(epsa->getResult()->lens[0]);
  auto Mabij(
    Tcc<TE>::template tensor<F>(std::vector<size_t>({Nv,Nv,No,No}),"Mabij")
  );

  COMPILE(
  //reconstruct mp2 amplitudes on-the-fly
    (*Mabij)["abij"] <<= map<F>(fromReal<F>, (*epsi)["i"]),
    (*Mabij)["abij"] +=  map<F>(fromReal<F>, (*epsi)["j"]),
    (*Mabij)["abij"] -=  map<F>(fromReal<F>, (*epsa)["a"]),
    (*Mabij)["abij"] -=  map<F>(fromReal<F>, (*epsa)["b"]),

    (*Mabij)["abij"] <<=
      map<F>(conj<F>, (*Vabij)["abij"]) *
      map<F>(inverse<F>, (*Mabij)["abij"]),
  // calculate mp2 pair energies in fno basis
    (*mp2PairEnergiesFno)["ij"] <<= ( 2.0) * (*Mabij)["abij"] * (*Vabij)["abij"],
    (*mp2PairEnergiesFno)["ij"]  += (-1.0) * (*Mabij)["abij"] * (*Vabij)["abji"],
    (*eMp2)[""] <<= (*mp2PairEnergiesFno)["ij"],
    (*eMp2Cbs)[""] <<= (*mp2PairEnergiesCbs)["ij"],

  //ccsd amplitudes
    (*Tabij)["abij"] <<= (*Tpphh)["abij"],
    (*Tabij)["abij"]  += (*Tph)["ai"] * (*Tph)["bj"],
  // reevalaute ccsd energy
    (*eCcsd)[""] <<= ( 2.0) * (*Tabij)["abij"] * (*Vabij)["abij"],
    (*eCcsd)[""]  += (-1.0) * (*Tabij)["abij"] * (*Vabij)["abji"],
  // evaluate nominator
    (*gijccd)["ij"] <<= (*Dabij)["abij"] * (*Tabij)["abij"],
    (*gijmp2)["ij"] <<= (*Dabij)["abij"] * (*Mabij)["abij"],
  // divide by <ij|δ|ij>
    (*gijccd)["ij"] <<= (*gijccd)["ij"] * map<F>(inverse<F>, (*nij)["ij"]),
    (*gijmp2)["ij"] <<= (*gijmp2)["ij"] * map<F>(inverse<F>, (*nij)["ij"]),

  // final multiplicative factor for the cbsEigenEnergies
    (*geff)["ij"] <<= map<Real<>>(real<F>, (*gijccd)["ij"]),
    (*geff)["ij"]  += map<Real<>>(real<F>, (*gijmp2)["ij"]),
    (*geff)["ij"]  += map<Real<>>(real<F>, (*gijccd)["ij"])
                    * map<Real<>>(real<F>, (*gijccd)["ij"]),

  // construct ΔEmp2 and scale with geff
    (*mp2PairEnergiesCbs)["ij"] +=
      (-1.0) * map<Real<>>(real<F>, (*mp2PairEnergiesFno)["ij"]),
    (*deltaEppl)[""] <<= (*geff)["ij"] * (*mp2PairEnergiesCbs)["ij"]

  )->execute();

  Real<> Emp2Cbs(eMp2Cbs->read());
  Real<> deltaPsPpl(deltaEppl->read());
  Real<> Emp2(real<F>(eMp2->read()));
  Real<> Eccsd(real<F>(eCcsd->read()));

  OUT() << "ccsd fno:         " << Eccsd << "\n";
  OUT() << "mp2 fno:          " << Emp2 << "\n";
  OUT() << "Δmp2:             " << Emp2Cbs - Emp2 << "\n";
  OUT() << "ps-ppl:           " << deltaPsPpl << "\n";
  OUT() << "corrected Energy: " << Eccsd - Emp2 + Emp2Cbs + deltaPsPpl << std::endl;

//  result->setValue<Ptr<Tr>>("geff", geff);
//  result->setValue<Ptr<T>>("gijccd", nij);
//  result->setValue<Ptr<Tr>>("mp2p", mp2PairEnergies);

  return true;
}
