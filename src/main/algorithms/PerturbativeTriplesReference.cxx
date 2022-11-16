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
#include <algorithms/CompleteRenormalized.hpp>
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
  if (Cc4s::dryRun) {
    using TE = DefaultDryTensorEngine;
    (
      result = run<Real<>,TE>(arguments)
    ) || (
      result = run<Complex<>,TE>(arguments)
    );
  } else {
    using TE = DefaultTensorEngine;
    (
      result = run<Real<>,TE>(arguments)
    ) || (
      result = run<Complex<>,TE>(arguments)
    );
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

  const bool completeRenormalized
     = arguments->getValue<int>("completeRenormalized", 0) == 1;

  auto eigenEnergies(
    arguments->getPtr<TensorSet<Real<>,TE>>("slicedEigenEnergies")
  );
  auto epsh(eigenEnergies->get("h"));
  auto epsp(eigenEnergies->get("p"));


  // piecuch three body amplitudes
  auto M( Tcc<TE>::template tensor<F>("M"));
  auto Z( Tcc<TE>::template tensor<F>("Z"));
  auto T( Tcc<TE>::template tensor<F>("T"));
  auto S( Tcc<TE>::template tensor<F>("S"));
  auto D( Tcc<TE>::template tensor<F>("D"));
  auto R( Tcc<TE>::template tensor<F>("R"));
  auto E( Tcc<TE>::template tensor<F>("E"));
  auto fromReal( [](Real<> x) { return F(x); } );
  auto inverse( [](F x) { return 1.0 / x; } );


  // deal with completeRenormalized
  if (completeRenormalized) {
    auto intermediates
       = cr::getCompleteRenormalized<F, TE>(coulombIntegrals, amplitudes);
    auto Jppph = intermediates->get("ppph");
    auto Jhphh = intermediates->get("hphh");
    COMPILE(
      (*M)["abcijk"] <<=          (*Jppph)["bcdk"] * (*Tpphh)["adij"],
      (*M)["abcijk"]  += (-1.0) * map<F>(conj<F>, (*Jhphh)["lcjk"])
                                * (*Tpphh)["abil"],
      (*Z)["abcijk"] <<= (*M)["abcijk"],
      (*Z)["abcijk"]  += (*M)["bacjik"],
      (*Z)["abcijk"]  += (*M)["acbikj"],
      (*Z)["abcijk"]  += (*M)["cbakji"],
      (*Z)["abcijk"]  += (*M)["cabkij"],
      (*Z)["abcijk"]  += (*M)["bcajki"],

      (*M)["abcijk"] <<= (*Z)["abcijk"]
    )->execute();
  }


  COMPILE(

    (*T)["abcijk"]  <<=          (*Vppph)["bcdk"] * (*Tpphh)["adij"],
    (*T)["abcijk"]   += (-1.0) * map<F>(conj<F>, (*Vhhhp)["jklc"]) * (*Tpphh)["abil"],
    (*Z)["abcijk"]  <<= (*T)["abcijk"],
    (*Z)["abcijk"]   += (*T)["bacjik"],
    (*Z)["abcijk"]   += (*T)["acbikj"],
    (*Z)["abcijk"]   += (*T)["cbakji"],
    (*Z)["abcijk"]   += (*T)["cabkij"],
    (*Z)["abcijk"]   += (*T)["bcajki"],

    (*T)["abcijk"]  <<= (*Z)["abcijk"],

    (*S)["abcijk"]  <<= (*Z)["abcijk"],
    (*S)["abcijk"]   += (*Tph)["ai"] * (*Vpphh)["bcjk"],
    (*S)["abcijk"]   += (*Tph)["bj"] * (*Vpphh)["acik"],
    (*S)["abcijk"]   += (*Tph)["ck"] * (*Vpphh)["abij"],

    (*Z)["abcijk"] <<= (4.0/3.0) * (*S)["abcijk"],
    (*Z)["abcijk"]  += (-2.0)    * (*S)["acbijk"],
    (*Z)["abcijk"]  += (2.0/3.0) * (*S)["bcaijk"],

    //just to tell tcc the dimensions of D
    (*D)["abcijk"] <<= (*T)["abcijk"],
    (*D)["abcijk"] <<= ( 1.0) * map<F>(fromReal, (*epsh)["i"]),
    (*D)["abcijk"]  += ( 1.0) * map<F>(fromReal, (*epsh)["j"]),
    (*D)["abcijk"]  += ( 1.0) * map<F>(fromReal, (*epsh)["k"]),
    (*D)["abcijk"]  += (-1.0) * map<F>(fromReal, (*epsp)["a"]),
    (*D)["abcijk"]  += (-1.0) * map<F>(fromReal, (*epsp)["b"]),
    (*D)["abcijk"]  += (-1.0) * map<F>(fromReal, (*epsp)["c"]),

    (*Z)["abcijk"]
      <<= map<F>(conj<F>, (*Z)["abcijk"]) * map<F>(inverse, (*D)["abcijk"]),

    (*T)["abcijk"]
      <<= map<F>(conj<F>, (*T)["abcijk"]) * map<F>(inverse, (*D)["abcijk"]),

    (*R)["abcijk"]  <<= (*Tph)["ai"] * (*Vpphh)["bcjk"],
    (*R)["abcijk"]   += (*Tph)["bj"] * (*Vpphh)["acik"],
    (*R)["abcijk"]   += (*Tph)["ck"] * (*Vpphh)["abij"],

    (*S)["abcijk"] <<= (4./3) * (*R)["abcijk"],
    (*S)["abcijk"]  += (-2.0) * (*R)["acbijk"],
    (*S)["abcijk"]  += (2./3) * (*R)["bcaijk"],

    (*S)["abcijk"]
      <<= map<F>(conj<F>, (*S)["abcijk"]) * map<F>(inverse, (*D)["abcijk"])




  )->execute();

  if (completeRenormalized) {
    COMPILE(
      (*E)[""] <<= (*M)["abcijk"] * (*Z)["abcijk"]
    )->execute();
  } else {
    COMPILE(
      (*E)[""] <<= (*T)["abcijk"] * (*Z)["abcijk"] * (*D)["abcijk"]
    )->execute();
  }

  /* DUPLICATION WARNING:
   * --------------------
   *
   * If you change something in the output of this algorithm
   * be sure to change accordingly the output of the
   * PerturbativeTriples, as they should be
   * symmetric.
   *
   */

  F eTriples(E->read());

  if (completeRenormalized) {
    double denominator = cr::getDenominator<F, TE>(amplitudes, T, S);
    OUT() << "CR-(T) correlation energy: "
          << std::setprecision(10) << std::setw(20)
          << real(eTriples) / denominator << std::endl;
    OUT() << "CR-(T) denominator: "
          << std::setprecision(10) << std::setw(27)
          << denominator << std::endl;
  } else {
    OUT() << "(T) correlation energy: "
          << std::setprecision(10) << std::setw(23)
          << real(eTriples) << std::endl;
  }

  auto energy(New<MapNode>(SOURCE_LOCATION));
  energy->setValue("correlation", real(eTriples));
  energy->setValue("unit", epsh->inspect()->getUnit());

  if (arguments->isGiven("mp2PairEnergies")) {
    const Real<>
      triples_star = computeCssdPtStar<F, TE>(arguments, real(eTriples));
    OUT() << "(T*)-Bsie energy correction: "
          << std::setprecision(15) << std::setw(23)
          << triples_star << std::endl;
    energy->setValue("starCorrection", real(triples_star));
  }

  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("energy") = energy;
  return result;
}

