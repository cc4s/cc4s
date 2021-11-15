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
  auto Vabij(coulombSlices->getValue<Ptr<TensorRecipe<F,TE>>>("pphh"));
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

  auto eigenEnergies(arguments->getMap("slicedEigenEnergies"));
  auto energySlices(eigenEnergies->getMap("slices"));
  auto epsi(energySlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("h"));
  auto epsa(energySlices->getValue<Ptr<TensorRecipe<Real<>,TE>>>("p"));

  auto No(epsi->getResult()->lens[0]);
  auto Nv(epsa->getResult()->lens[0]);
  auto Dabij(
    Tcc<TE>::template tensor<F>(std::vector<size_t>({Nv,Nv,No,No}),"Dabij")
  );
  auto direct( Tcc<TE>::template tensor<F>("D") );
  auto exchange( Tcc<TE>::template tensor<F>("X") );
  auto toEigenUnits = pow(eigenEnergies->getValue<Real<>>("unit"),2.0) /
    pow(coulombIntegrals->getValue<Real<>>("unit"),2.0);
  OUT() << "Contracting MP2 energy..." << std::endl;
  COMPILE(
    (*Dabij)["abij"] <<= map<F>([](Real<> eps) {return F(eps);}, (*epsa)["a"]),
    (*Dabij)["abij"] +=  map<F>([](Real<> eps) {return F(eps);}, (*epsa)["b"]),
    (*Dabij)["abij"] -=  map<F>([](Real<> eps) {return F(eps);}, (*epsi)["i"]),
    (*Dabij)["abij"] -=  map<F>([](Real<> eps) {return F(eps);}, (*epsi)["j"]),
    (*Dabij)["abij"] <<=
      map<F>(conj<F>, (*Vabij)["abij"]) *
      map<F>([](F delta) { return F(1)/delta; }, (*Dabij)["abij"]),
    (*direct)[""] <<=
      toEigenUnits * -0.5*spins*spins * (*Vabij)["abij"] * (*Dabij)["abij"],
    (*exchange)[""] <<=
      toEigenUnits * +0.5*spins * (*Vabij)["abji"] * (*Dabij)["abij"]
  )->execute();

  F D(direct->read()), X(exchange->read());
  OUT() << "energy= " << D+X << std::endl;
  OUT() << "\tdirect= " << D << std::endl;
  OUT() << "\texchange= " << X << std::endl;
  
  auto energy(New<MapNode>(SOURCE_LOCATION));
  energy->setValue<Real<>>("direct", real<F>(D));
  energy->setValue<Real<>>("exchange", real<F>(X));
  energy->setValue<Real<>>("value", real<F>(D+X));
  energy->setValue<Real<>>("unit", eigenEnergies->getValue<Real<>>("unit"));
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("energy") = energy;
  return result;
}

