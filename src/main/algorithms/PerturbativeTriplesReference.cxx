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

  auto amplitudes = arguments->getPtr<TensorUnion<F,TE>>("amplitudes");
  auto Tph( amplitudes->get(0) );
  auto Tpphh( amplitudes->get(1) );

  auto coulombIntegrals(arguments->getMap("coulombIntegrals"));
  auto coulombSlices(coulombIntegrals->getMap("slices"));
  auto Vpphh(coulombSlices->getPtr<TensorExpression<F,TE>>("pphh"));
  auto Vppph(coulombSlices->getPtr<TensorExpression<F,TE>>("ppph"));
  auto Vhhhp(coulombSlices->getPtr<TensorExpression<F,TE>>("hhhp"));

  auto eigenEnergies(arguments->getMap("slicedEigenEnergies"));
  auto energySlices(eigenEnergies->getMap("slices"));
  auto epsi(energySlices->getPtr<TensorExpression<Real<>,TE>>("h"));
  auto epsa(energySlices->getPtr<TensorExpression<Real<>,TE>>("p"));


  auto Z( Tcc<TE>::template tensor<F>("Z"));
  auto T( Tcc<TE>::template tensor<F>("T"));
  auto S( Tcc<TE>::template tensor<F>("S"));
  auto E( Tcc<TE>::template tensor<F>("E"));
  auto fromReal( [](Real<> x) { return F(x); } );
  auto inverse( [](F x) { return 1.0 / x; } );
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

    (*S)["abcijk"] <<= ( 1.0) * map<F>(fromReal, (*epsi)["i"]),
    (*S)["abcijk"]  += ( 1.0) * map<F>(fromReal, (*epsi)["j"]),
    (*S)["abcijk"]  += ( 1.0) * map<F>(fromReal, (*epsi)["k"]),
    (*S)["abcijk"]  += (-1.0) * map<F>(fromReal, (*epsa)["a"]),
    (*S)["abcijk"]  += (-1.0) * map<F>(fromReal, (*epsa)["b"]),
    (*S)["abcijk"]  += (-1.0) * map<F>(fromReal, (*epsa)["c"]),

    (*Z)["abcijk"] <<= (*Z)["abcijk"] * map<F>(inverse, (*S)["abcijk"]),
    (*E)[""] <<= (*T)["abcijk"] * (*Z)["abcijk"],
    (*E)[""]  += (*T)["bacjik"] * (*Z)["abcijk"],
    (*E)[""]  += (*T)["acbikj"] * (*Z)["abcijk"],
    (*E)[""]  += (*T)["cbakji"] * (*Z)["abcijk"],
    (*E)[""]  += (*T)["cabkij"] * (*Z)["abcijk"],
    (*E)[""]  += (*T)["bcajki"] * (*Z)["abcijk"]
  )->execute();

  F eTriples(E->read());
  OUT() << "(T) correlation energy: "
        << std::setprecision(15) << std::setw(23)
        << real(eTriples) << std::endl;

  auto energy(New<MapNode>(SOURCE_LOCATION));
  energy->setValue("triples", real(eTriples));
  energy->setValue("unit", eigenEnergies->getValue<Real<>>("unit"));
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("energy") = energy;
  return result;
}

