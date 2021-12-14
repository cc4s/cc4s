/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <algorithms/SecondOrderPerturbationTheory.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(SecondOrderPerturbationTheory)

Ptr<MapNode> SecondOrderPerturbationTheory::run(const Ptr<MapNode> &arguments) {
  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto scalarType(coulombIntegrals->getValue<std::string>("scalarType"));
  // multiplex calls to template methods
  if (Cc4s::options->dryRun) {
    using TE = DefaultDryTensorEngine;
    if (scalarType == TypeTraits<Real<>>::getName()) {
      return calculateMp2Energy<Real<>,TE>(arguments);
    } else if (scalarType == TypeTraits<Complex<>>::getName()) {
      return calculateMp2Energy<Complex<>,TE>(arguments);
    }
  } else {
    using TE = DefaultTensorEngine;
    if (scalarType == TypeTraits<Real<>>::getName()) {
      return calculateMp2Energy<Real<>,TE>(arguments);
    } else if (scalarType == TypeTraits<Complex<>>::getName()) {
      return calculateMp2Energy<Complex<>,TE>(arguments);
    }
  }
  ASSERT_LOCATION(
    false, "unsupported orbitals type '" + scalarType + "'",
    coulombIntegrals->get("scalarType")->sourceLocation
  );
}

template <typename F, typename TE>
Ptr<MapNode> SecondOrderPerturbationTheory::calculateMp2Energy(
  const Ptr<MapNode> &arguments
) {
  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto coulombSlices(coulombIntegrals->getMap("slices"));
  auto Vpphh(coulombSlices->getPtr<TensorExpression<F,TE>>("pphh"));
  auto orbitalType(
    coulombIntegrals->getMap(
      "indices"
    )->getMap("orbital")->getValue<std::string>("type")
  );
  Real<> spins;
  if (orbitalType == "spatial") {
    spins = 2;
  } else if (orbitalType == "spin") {
    spins = 1;
  } else {
    ASSERT_LOCATION(
      false, "unsupported orbital type '" + orbitalType + "'",
      coulombIntegrals->getMap(
        "indices"
      )->getMap("orbital")->get("type")->sourceLocation
    );
  }

  auto fockOperator(getFockOperator<F,TE>(arguments));
  auto fockSlices(fockOperator->getMap("slices"));
  auto fph(fockSlices->template getPtr<TensorExpression<F,TE>>("ph"));

  auto eigenEnergies(arguments->getMap("slicedEigenEnergies"));
  auto energySlices(eigenEnergies->getMap("slices"));
  auto epsi(energySlices->getPtr<TensorExpression<Real<>,TE>>("h"));
  auto epsa(energySlices->getPtr<TensorExpression<Real<>,TE>>("p"));

  auto No(epsi->inspect()->getLen(0));
  auto Nv(epsa->inspect()->getLen(0));
  auto Dph(
    Tcc<TE>::template tensor<F>(std::vector<size_t>({Nv,No}),"Dph")
  );
  auto Dpphh(
    Tcc<TE>::template tensor<F>(std::vector<size_t>({Nv,Nv,No,No}),"Dpphh")
  );
  auto singles( Tcc<TE>::template tensor<F>("S") );
  auto direct( Tcc<TE>::template tensor<F>("D") );
  auto exchange( Tcc<TE>::template tensor<F>("X") );
  auto toEigenUnits(
    pow(eigenEnergies->getValue<Real<>>("unit"),2.0)
    / pow(coulombIntegrals->getValue<Real<>>("unit"),2.0)
  );
  OUT() << "Contracting MP2 energy..." << std::endl;
  COMPILE(
    (*Dph)["ai"] <<= map<F>([](Real<> eps) {return F(eps);}, (*epsa)["a"]),
    (*Dph)["ai"] -=  map<F>([](Real<> eps) {return F(eps);}, (*epsi)["i"]),
    (*Dpphh)["abij"] <<= (*Dph)["ai"],
    (*Dpphh)["abij"] +=  (*Dph)["bj"],
    // first-order singles amplitudes
    (*Dph)["ai"] <<=
      map<F>(conj<F>, (*fph)["ai"]) *
      map<F>([](F delta) { return F(1/real(delta)); }, (*Dph)["ai"]),
    // first-order doubles amplitudes
    (*Dpphh)["abij"] <<=
      map<F>(conj<F>, (*Vpphh)["abij"]) *
      map<F>([](F delta) { return F(1/real(delta)); }, (*Dpphh)["abij"]),
    (*singles)[""] <<=
      toEigenUnits * -spins * (*fph)["ai"] * (*Dph)["ai"],
    (*direct)[""] <<=
      toEigenUnits * -0.5*spins*spins * (*Vpphh)["abij"] * (*Dpphh)["abij"],
    (*exchange)[""] <<=
      toEigenUnits * +0.5*spins * (*Vpphh)["abji"] * (*Dpphh)["abij"]
  )->execute();

  F S(singles->read()), D(direct->read()), X(exchange->read());
  OUT() << "correlation energy: " << S+D+X << std::endl;
  OUT() << "  singles:  " << S << std::endl;
  OUT() << "  direct:   " << D << std::endl;
  OUT() << "  exchange: " << X << std::endl;

  auto energy(New<MapNode>(SOURCE_LOCATION));
  energy->setValue("singles", real(S));
  energy->setValue("direct", real(D));
  energy->setValue("exchange", real(X));
  energy->setValue("correlation", real(S+D+X));
  energy->setValue("secondOrderSingles", real(S));
  energy->setValue("secondOrderDirect", real(D));
  energy->setValue("secondOrderExchange", real(X));
  energy->setValue("secondOrder", real(S+D+X));
  energy->setValue("unit", eigenEnergies->getValue<Real<>>("unit"));
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("energy") = energy;
  return result;
}

template <typename F, typename TE>
Ptr<MapNode> SecondOrderPerturbationTheory::getFockOperator(
  const Ptr<MapNode> &arguments
) {
  if (arguments->isGiven("slicedFockOperator")) {
    return arguments->getMap("slicedFockOperator");
  }
  // else assume canonical calculation and compute from eigenenergies
  auto eigenEnergies(arguments->getMap("slicedEigenEnergies"));
  auto energySlices(eigenEnergies->getMap("slices"));
  auto epsi(energySlices->getPtr<TensorExpression<Real<>,TE>>("h"));
  auto epsa(energySlices->getPtr<TensorExpression<Real<>,TE>>("p"));
  auto No(epsi->inspect()->getLen(0));
  auto Nv(epsa->inspect()->getLen(0));
  
  auto fhh(Tcc<TE>::template tensor<F>(std::vector<Natural<>>({No,No}),"fhh"));
  auto fph(Tcc<TE>::template tensor<F>(std::vector<Natural<>>({Nv,No}),"fph"));
  auto fhp(Tcc<TE>::template tensor<F>(std::vector<Natural<>>({No,Nv}),"fhp"));
  auto fpp(Tcc<TE>::template tensor<F>(std::vector<Natural<>>({Nv,Nv}),"fpp"));

  // off-diagonal slices are zero tensors
  // diagonal slices have eigenenergies on diagonal
  COMPILE(
    (*fhh)["ii"] <<= map<F>([](Real<> eps) {return F(eps);}, (*epsi)["i"]),
    (*fpp)["aa"] <<= map<F>([](Real<> eps) {return F(eps);}, (*epsa)["a"])
  )->execute();

  auto result( New<MapNode>(SOURCE_LOCATION) );
  result->setValue("unit", eigenEnergies->getValue<Real<>>("unit"));
  result->setValue("scalarType", TypeTraits<F>::getName());
  auto slices( New<MapNode>(SOURCE_LOCATION) );
  slices->setPtr("hh", fhh);
  slices->setPtr("ph", fph);
  slices->setPtr("hp", fhp);
  slices->setPtr("pp", fpp);
  result->get("slices") = slices;
  return result;
}

