#include <algorithms/PerturbativeTriplesReference.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <math/TensorUnion.hpp>


using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(PerturbativeTriplesReference)

Ptr<MapNode> PerturbativeTriplesReference::run(const Ptr<MapNode> &arguments) {
  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto scalarType(coulombIntegrals->getValue<std::string>("scalarType"));
  // multiplex calls to template methods
  if (Cc4s::options->dryRun) {
    using TE = DefaultDryTensorEngine;
    if (scalarType == TypeTraits<Real<>>::getName()) {
      return calculateTriplesEnergy<Real<>,TE>(arguments);
    } else if (scalarType == TypeTraits<Complex<>>::getName()) {
      return calculateTriplesEnergy<Complex<>,TE>(arguments);
    }
  } else {
    using TE = DefaultTensorEngine;
    if (scalarType == TypeTraits<Real<>>::getName()) {
      return calculateTriplesEnergy<Real<>,TE>(arguments);
    } else if (scalarType == TypeTraits<Complex<>>::getName()) {
      return calculateTriplesEnergy<Complex<>,TE>(arguments);
    }
  }
  ASSERT_LOCATION(
    false, "unsupported orbitals type '" + scalarType + "'",
    coulombIntegrals->get("scalarType")->sourceLocation
  );
}

template <typename F, typename TE>
Ptr<MapNode> PerturbativeTriplesReference::calculateTriplesEnergy(
  const Ptr<MapNode> &arguments
) {

  auto amplitudes = arguments->getValue<Ptr<const TensorUnion<F,TE>>>("amplitudes");
  auto Tph( amplitudes->get(0) );
  auto Tpphh( amplitudes->get(1) );

  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto coulombSlices(coulombIntegrals->getMap("slices"));
  auto Vpphh(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("pphh"));
  auto Vppph(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("ppph"));
  auto Vhhhp(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("hhhp"));

  auto eigenEnergies(arguments->getMap("slicedEigenEnergies"));
  auto energySlices(eigenEnergies->getMap("slices"));
  auto epsi(energySlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("h"));
  auto epsa(energySlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("p"));


  auto Z( Tcc<TE>::template tensor<F>("Z"));
  auto T( Tcc<TE>::template tensor<F>("T"));
  auto S( Tcc<TE>::template tensor<F>("S"));
  auto E( Tcc<TE>::template tensor<F>("E"));
  COMPILE(
    (*T)["abcijk"]  <<=          (*Vppph)["bcdk"] * (*Tpphh)["adij"],
    (*T)["abcijk"]   += (-1.0) * (*Vhhhp)["jklc"] * (*Tpphh)["abil"],
    (*S)["abcijk"]  <<= (*T)["abcijk"],
    (*S)["abcijk"]   += ( 0.5) * (*Tph)["ai"] * (*Vpphh)["bcjk"],

    (*Z)["abcijk"] <<= (+8.0) * (*S)["abcijk"],
    (*Z)["abcijk"]  += (-4.0) * (*S)["acbijk"],
    (*Z)["abcijk"]  += (-4.0) * (*S)["bacijk"],
    (*Z)["abcijk"]  += (+2.0) * (*S)["bcaijk"],
    (*Z)["abcijk"]  += (+2.0) * (*S)["cabijk"],
    (*Z)["abcijk"]  += (-4.0) * (*S)["cbaijk"],

    (*S)["abcijk"] <<= ( 1.0) * map<F>(fromReal<F>, (*epsi)["i"]),
    (*S)["abcijk"]  += ( 1.0) * map<F>(fromReal<F>, (*epsi)["j"]),
    (*S)["abcijk"]  += ( 1.0) * map<F>(fromReal<F>, (*epsi)["k"]),
    (*S)["abcijk"]  += (-1.0) * map<F>(fromReal<F>, (*epsa)["a"]),
    (*S)["abcijk"]  += (-1.0) * map<F>(fromReal<F>, (*epsa)["b"]),
    (*S)["abcijk"]  += (-1.0) * map<F>(fromReal<F>, (*epsa)["c"]),

    (*Z)["abcijk"] <<= (*Z)["abcijk"] * map<F>(inverse<F>, (*S)["abcijk"]),
    (*E)[""] <<= (*T)["abcijk"] * (*Z)["abcijk"],
    (*E)[""]  += (*T)["bacjik"] * (*Z)["abcijk"],
    (*E)[""]  += (*T)["acbikj"] * (*Z)["abcijk"],
    (*E)[""]  += (*T)["cbakji"] * (*Z)["abcijk"],
    (*E)[""]  += (*T)["cabkij"] * (*Z)["abcijk"],
    (*E)[""]  += (*T)["bcajki"] * (*Z)["abcijk"]
  )->execute();

  F eTriples(E->read());
  OUT() << "Perturbative Triples Energy: " << real<F>(eTriples) << std::endl;

  auto energy(New<MapNode>(SOURCE_LOCATION));
  energy->setValue<Real<>>("triples", real<F>(eTriples));
  energy->setValue<Real<>>("unit", eigenEnergies->getValue<Real<>>("unit"));
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("energy") = energy;
  return result;
}

