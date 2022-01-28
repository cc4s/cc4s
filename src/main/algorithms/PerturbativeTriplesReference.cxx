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

#include <algorithms/PerturbativeTriplesStar.hpp>
#include <algorithms/PerturbativeTriplesReference.hpp>
#include <tcc/Tcc.hpp>
#include <TensorSet.hpp>
#include <MathFunctions.hpp>
#include <Log.hpp>
#include <Exception.hpp>
#include <Cc4s.hpp>


using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(PerturbativeTriplesReference)

Ptr<MapNode> PerturbativeTriplesReference::run(const Ptr<MapNode> &arguments) {
  Ptr<MapNode> result;
  // multiplex calls to template methods
  if (Cc4s::options->dryRun) {
    using TE = DefaultDryTensorEngine;
    (result = run<Real<>,TE>(arguments))
      || (result = run<Complex<>,TE>(arguments));
  } else {
    using TE = DefaultTensorEngine;
    (result = run<Real<>,TE>(arguments))
      || (result = run<Complex<>,TE>(arguments));
  }
  ASSERT_LOCATION(
    result, "unsupported tensor type as 'amplitudes'",
    arguments->sourceLocation
  );
  return result;
}

template <typename F, typename TE>
Ptr<MapNode> PerturbativeTriplesReference::run(const Ptr<MapNode> &arguments) {
  auto amplitudes = arguments->getPtr<TensorSet<F,TE>>("amplitudes");
  if (!amplitudes) return nullptr;
  auto Tph( amplitudes->get("ph") );
  auto Tpphh( amplitudes->get("pphh") );

  auto coulombIntegrals(arguments->getPtr<TensorSet<F,TE>>("coulombIntegrals"));
  auto Vpphh(coulombIntegrals->get("pphh"));
  auto Vppph(coulombIntegrals->get("ppph"));
  auto Vhhhp(coulombIntegrals->get("hhhp"));

  auto eigenEnergies(
    arguments->getPtr<TensorSet<Real<>,TE>>("slicedEigenEnergies")
  );
  auto epsh(eigenEnergies->get("h"));
  auto epsp(eigenEnergies->get("p"));


  auto Z( Tcc<TE>::template tensor<F>("Z"));
  auto T( Tcc<TE>::template tensor<F>("T"));
  auto S( Tcc<TE>::template tensor<F>("S"));
  auto E( Tcc<TE>::template tensor<F>("E"));
  auto fromReal( [](Real<> x) { return F(x); } );
  auto inverse( [](F x) { return 1.0 / x; } );
  COMPILE(
    (*T)["abcijk"]  <<=          (*Vppph)["bcdk"] * (*Tpphh)["adij"],
    (*T)["abcijk"]   += (-1.0) * map<F>(conj<F>, (*Vhhhp)["jklc"]) * (*Tpphh)["abil"],
    (*S)["abcijk"]  <<= (*T)["abcijk"],
    (*S)["abcijk"]   += ( 0.5) * (*Tph)["ai"] * (*Vpphh)["bcjk"],

    (*Z)["abcijk"] <<= (+8.0) * (*S)["abcijk"],
    (*Z)["abcijk"]  += (-4.0) * (*S)["acbijk"],
    (*Z)["abcijk"]  += (-4.0) * (*S)["bacijk"],
    (*Z)["abcijk"]  += (+2.0) * (*S)["bcaijk"],
    (*Z)["abcijk"]  += (+2.0) * (*S)["cabijk"],
    (*Z)["abcijk"]  += (-4.0) * (*S)["cbaijk"],

    (*S)["abcijk"] <<= ( 1.0) * map<F>(fromReal, (*epsh)["i"]),
    (*S)["abcijk"]  += ( 1.0) * map<F>(fromReal, (*epsh)["j"]),
    (*S)["abcijk"]  += ( 1.0) * map<F>(fromReal, (*epsh)["k"]),
    (*S)["abcijk"]  += (-1.0) * map<F>(fromReal, (*epsp)["a"]),
    (*S)["abcijk"]  += (-1.0) * map<F>(fromReal, (*epsp)["b"]),
    (*S)["abcijk"]  += (-1.0) * map<F>(fromReal, (*epsp)["c"]),

    (*Z)["abcijk"] <<= (*Z)["abcijk"] * map<F>(inverse, (*S)["abcijk"]),
    (*E)[""] <<= map<F>(conj<F>, (*T)["abcijk"]) * (*Z)["abcijk"],
    (*E)[""]  += map<F>(conj<F>, (*T)["bacjik"]) * (*Z)["abcijk"],
    (*E)[""]  += map<F>(conj<F>, (*T)["acbikj"]) * (*Z)["abcijk"],
    (*E)[""]  += map<F>(conj<F>, (*T)["cbakji"]) * (*Z)["abcijk"],
    (*E)[""]  += map<F>(conj<F>, (*T)["cabkij"]) * (*Z)["abcijk"],
    (*E)[""]  += map<F>(conj<F>, (*T)["bcajki"]) * (*Z)["abcijk"]
  )->execute();

  F eTriples(E->read());
  OUT() << "(T) correlation energy: "
        << std::setprecision(15) << std::setw(23)
        << real(eTriples) << std::endl;

  auto energy(New<MapNode>(SOURCE_LOCATION));
  energy->setValue("correlation", real(eTriples));
  energy->setValue("unit", epsh->inspect()->getUnit());

  if (arguments->isGiven("mp2PairEnergies")) {
    const Real<>
      triples_star = computeCssdPtStar<F, TE>(arguments, real(eTriples));
    OUT() << "(T*) correlation energy: "
          << std::setprecision(15) << std::setw(23)
          << triples_star << std::endl;
    energy->setValue("starCorrelation", real(triples_star));
  }

  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("energy") = energy;
  return result;
}

