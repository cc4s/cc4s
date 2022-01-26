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
#include <tcc/Tcc.hpp>
#include <MathFunctions.hpp>
#include <Log.hpp>
#include <Exception.hpp>
#include <Cc4s.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(SecondOrderPerturbationTheory)

Ptr<MapNode> SecondOrderPerturbationTheory::run(const Ptr<MapNode> &arguments) {
  // multiplex calls to template methods
  Ptr<MapNode> result;
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
    result, "unsupported tensor type as 'operator'",
    arguments->sourceLocation
  );
  return result;
}

template <typename F, typename TE>
Ptr<MapNode> SecondOrderPerturbationTheory::run(
  const Ptr<MapNode> &arguments
) {
  auto coulombIntegrals(arguments->getPtr<TensorSet<F,TE>>("coulombIntegrals"));
  if (!coulombIntegrals) return nullptr;
  auto Vpphh(coulombIntegrals->get("pphh"));
  // TODO: find number of spin properties
  Real<> degeneracy(2);

  auto fockOperator(getFockOperator<F,TE>(arguments));
  auto fph(fockOperator->get("ph"));

  auto eigenEnergies(
    arguments->getPtr<TensorSet<Real<>,TE>>("slicedEigenEnergies")
  );
  auto epsh(eigenEnergies->get("h"));
  auto epsp(eigenEnergies->get("p"));

  auto No(epsh->inspect()->getLen(0));
  auto Nv(epsp->inspect()->getLen(0));
  auto Dph(
    Tcc<TE>::template tensor<F>(std::vector<Natural<>>({Nv,No}),"Dph")
  );
  auto Dpphh(
    Tcc<TE>::template tensor<F>(std::vector<Natural<>>({Nv,Nv,No,No}),"Dpphh")
  );
  auto singles( Tcc<TE>::template tensor<F>("S") );
  auto direct( Tcc<TE>::template tensor<F>("D") );
  auto exchange( Tcc<TE>::template tensor<F>("X") );
  OUT() << "Contracting second order energy..." << std::endl;
  COMPILE(
    (*Dph)["ai"] <<= map<F>([](Real<> eps) {return F(eps);}, (*epsp)["a"]),
    (*Dph)["ai"] -=  map<F>([](Real<> eps) {return F(eps);}, (*epsh)["i"]),
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
      -degeneracy * (*fph)["ai"] * (*Dph)["ai"],
    (*direct)[""] <<=
      -0.5*degeneracy*degeneracy * (*Vpphh)["abij"] * (*Dpphh)["abij"],
    (*exchange)[""] <<=
      +0.5*degeneracy * (*Vpphh)["abji"] * (*Dpphh)["abij"]
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
  // TODO: use unit of scalar as determined by tcc
  energy->setValue("unit", epsh->inspect()->getUnit());
  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("energy") = energy;
  return result;
}

template <typename F, typename TE>
Ptr<TensorSet<F,TE>> SecondOrderPerturbationTheory::getFockOperator(
  const Ptr<MapNode> &arguments
) {
  if (arguments->isGiven("slicedFockOperator")) {
    return arguments->getPtr<TensorSet<F,TE>>("slicedFockOperator");
  }
  // else assume canonical calculation and compute from eigenenergies
  auto eigenEnergies(
    arguments->getPtr<TensorSet<Real<>,TE>>("slicedEigenEnergies")
  );
  auto epsh(eigenEnergies->get("h"));
  auto epsp(eigenEnergies->get("p"));
  auto No(epsh->inspect()->getLen(0));
  auto Nv(epsp->inspect()->getLen(0));
  
  auto fhh(Tcc<TE>::template tensor<F>(std::vector<Natural<>>({No,No}),"fhh"));
  auto fph(Tcc<TE>::template tensor<F>(std::vector<Natural<>>({Nv,No}),"fph"));
  auto fhp(Tcc<TE>::template tensor<F>(std::vector<Natural<>>({No,Nv}),"fhp"));
  auto fpp(Tcc<TE>::template tensor<F>(std::vector<Natural<>>({Nv,Nv}),"fpp"));

  // off-diagonal slices are zero tensors
  // diagonal slices have eigenenergies on diagonal
  COMPILE(
    (*fhh)["ii"] <<= map<F>([](Real<> eps) {return F(eps);}, (*epsh)["i"]),
    (*fpp)["aa"] <<= map<F>([](Real<> eps) {return F(eps);}, (*epsp)["a"])
  )->execute();

  auto result( New<TensorSet<F,TE>>() );
  result->get("hh") = fhh;
  result->get("ph") = fph;
  result->get("hp") = fhp;
  result->get("pp") = fpp;
  return result;
}

