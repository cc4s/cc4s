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
  using T = Tensor<F, TE>;
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
  auto Tabij( Tcc<TE>::template tensor<F>("Tabij"));

  auto gijccd( Tcc<TE>::template tensor<F>("gijccd"));
  auto gijmp2( Tcc<TE>::template tensor<F>("gijmp2"));

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

  //ccsd amplitudes
    (*Tabij)["abij"] <<= (*Tpphh)["abij"],
    (*Tabij)["abij"]  += (*Tph)["ai"] * (*Tph)["bj"],
  // evaluate nominator
    (*gijccd)["ij"] <<= (*Dabij)["abij"] * (*Tabij)["abij"],
    (*gijmp2)["ij"] <<= (*Dabij)["abij"] * (*Mabij)["abij"],
  // divide by <ij|δ|ij>
    (*gijccd)["ij"] <<= (*gijccd)["ij"] * map<F>(inverse<F>, (*nij)["ij"]),
    (*gijmp2)["ij"] <<= (*gijmp2)["ij"] * map<F>(inverse<F>, (*nij)["ij"])

  )->execute();


  result->setValue<Ptr<Tensor<F,TE>>>("gijccd", gijccd);
  result->setValue<Ptr<Tensor<F,TE>>>("gijmp2", gijmp2);
  return true;
}
