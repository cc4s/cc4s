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

#include <algorithms/coupledcluster/Drccd.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/SharedPointer.hpp>
#include <util/Exception.hpp>

using namespace cc4s;

template <typename F, typename TE>
CoupledClusterMethodRegistrar<
  F,TE,Drccd<F,TE>
> Drccd<F,TE>::registrar_("Drccd");


template <typename F, typename TE>
Ptr<TensorSet<F,TE>> Drccd<F,TE>::getResiduum(
  const Ptr<TensorSet<F,TE>> &amplitudes
) {
  // read all required integrals
  auto coulombIntegrals(this->arguments->getMap("coulombIntegrals"));
  auto coulombSlices(coulombIntegrals->getMap("slices"));
  auto Vpphh(coulombSlices->template getPtr<TensorExpression<F,TE>>("pphh"));
  auto Vphhp(coulombSlices->template getPtr<TensorExpression<F,TE>>("phhp"));
  auto Vhhpp(coulombSlices->template getPtr<TensorExpression<F,TE>>("hhpp"));

  // get spins
  auto orbitalType(
    coulombIntegrals->getMap(
      "indices"
    )->getMap("orbital")->template getValue<std::string>("type")
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

  // construct residuum. Shape will be assumed upon first use.
  auto Rph( Tcc<TE>::template tensor<F>("Rph") );
  auto Rpphh( Tcc<TE>::template tensor<F>("Rpphh") );
  auto residuum(
    New<TensorSet<F,TE>>(
      std::vector<Ptr<TensorExpression<F,TE>>>({Rph, Rpphh}),
      std::vector<std::string>({"ai", "abij"})
    )
  );

  bool linearized(
    this->arguments->template getValue<bool>("linearized", false)
  );
  bool adjacentPairsExchange(
    this->arguments->template getValue<bool>("adjacentPairsExchange", false)
  );

  if (!amplitudes) {
    // no previous amplitudes given
    COMPILE(
      (*Rph)["ai"] <<= F(0.0) * (*Vpphh)["aaii"],
      (*Rpphh)["abij"] <<= (*Vpphh)["abij"]
    )->execute();
  } else {
    // TODO: check if given amplitudes contain expected parts
    // get amplitude parts
    auto Tph( amplitudes->get(0) );
    auto Tpphh( amplitudes->get(1) );
    Tph->inspect()->setName("Tph"); Tpphh->inspect()->setName("Tpphh");

    auto Whhpp( Tcc<TE>::template tensor<F>("Whhpp") );
    // for the remaining iterations compute the drCCD residuum
    COMPILE(
      (*Rph)["ai"] <<= F(0.0) * (*Vpphh)["aaii"],
      (*Rpphh)["abij"] <<= spins * (*Vphhp)["akic"] * (*Tpphh)["cbkj"],
      (*Rpphh)["abij"] += (*Rpphh)["baji"],
      (*Rpphh)["abij"] += (*Vpphh)["abij"],
      (linearized) ? (
        // linearized: nothing more to do
        Tcc<TE>::sequence()
      ) : (
        // otherwise: do quadratic contribution
        (*Whhpp)["ijab"] <<= spins * (*Vhhpp)["ijab"],
        (adjacentPairsExchange) ? (
          // adjacent pairs correction: also exchange holes in Whhpp
          (*Whhpp)["ijab"] -= (*Vhhpp)["jiab"],
          Tcc<TE>::sequence()
        ) : (
          // otherwise: do nothing else with Whhpp
          Tcc<TE>::sequence()
        ),
        // compute quadratic contribution
        (*Rpphh)["abij"] +=
          spins * (*Whhpp)["klcd"] * (*Tpphh)["acik"] * (*Tpphh)["dblj"],
        Tcc<TE>::sequence()
      )
    )->execute();
  }

  return residuum;
}

// instantiate
template class cc4s::Drccd<Real<64>, DefaultDryTensorEngine>;
template class cc4s::Drccd<Complex<64>, DefaultDryTensorEngine>;
template class cc4s::Drccd<Real<64>, DefaultTensorEngine>;
template class cc4s::Drccd<Complex<64>, DefaultTensorEngine>;

