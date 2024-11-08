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

#include <algorithms/CcsdChannels.hpp>
#include <MathFunctions.hpp>
#include <Log.hpp>
#include <SharedPointer.hpp>
#include <Exception.hpp>
#include <TensorSet.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(CcsdChannels)

#define _FORMAT(_fmt, ...)                                                     \
  ([&](void) -> std::string {                                                  \
    int _sz = std::snprintf(nullptr, 0, _fmt, __VA_ARGS__);                    \
    std::vector<char> _out(_sz + 1);                                           \
    std::snprintf(&_out[0], _out.size(), _fmt, __VA_ARGS__);                   \
    return std::string(_out.data());                                           \
  })()


Ptr<MapNode> CcsdChannels::run(const Ptr<MapNode> &arguments) {
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
    result, "Unsupported Amplitudes",
    arguments->sourceLocation
  );
  return result;
}


template <typename F, typename TE>
Ptr<MapNode> CcsdChannels::run(const Ptr<MapNode> &arguments) {

  auto Rph( Tcc<TE>::template tensor<F>("Rph") );
  auto Rpphh( Tcc<TE>::template tensor<F>("Rpphh") );

  auto amplitudes(
    arguments->template getPtr<TensorSet<F,TE>>("amplitudes")
  );

  if (!amplitudes) return nullptr;

  auto coulombIntegrals(
    arguments->template getPtr<TensorSet<F,TE>>("coulombIntegrals")
  );
  ASSERT_LOCATION(
    coulombIntegrals, "Missing Integrals!",
    arguments->sourceLocation
  );
  auto Vpphh(coulombIntegrals->get("pphh"));


  auto eigenEnergies(
    arguments->template getPtr<TensorSet<Real<>,TE>>("slicedEigenEnergies")
  );
  ASSERT_LOCATION(
    eigenEnergies, "Missing EigenEnergies!",
    arguments->sourceLocation
  );

  auto Tph( amplitudes->get("ph") );
  auto Tpphh( amplitudes->get("pphh") );
  Tph->inspect()->setName("Tph"); Tpphh->inspect()->setName("Tpphh");


  auto epsh(eigenEnergies->get("h"));
  auto epsp(eigenEnergies->get("p"));
  auto Fepsh(Tcc<TE>::template tensor<F>(epsh->inspect()->getLens(), "Fepsh"));
  auto Fepsp(Tcc<TE>::template tensor<F>(epsp->inspect()->getLens(), "Fepsp"));
  auto fromReal( [](Real<> eps) {return F(eps);} );
  auto Dph( Tcc<TE>::template tensor<F>("Dph") );
  auto Dpphh( Tcc<TE>::template tensor<F>("Dpphh") );
  COMPILE(
    (*Dph)["ai"] <<= (0.0) * (*Vpphh)["aaii"],
    (*Dpphh)["abij"] <<= (0.0) * (*Vpphh)["abij"],
    (*Fepsp)["a"] <<= map<F>(fromReal, (*epsp)["a"]),
    (*Fepsh)["i"] <<= map<F>(fromReal, (*epsh)["i"]),
    (*Dpphh)["abij"]  += (*Fepsp)["a"],
    (*Dpphh)["abij"]  += (*Fepsp)["b"],
    (*Dpphh)["abij"]  -= (*Fepsh)["i"],
    (*Dpphh)["abij"]  -= (*Fepsh)["j"],
    (*Dph)["ai"]  += (*Fepsp)["a"],
    (*Dph)["ai"]  -= (*Fepsh)["i"]
  )->execute();
  Dph->sum(F(1.0), "ai", "ai", [] (F delta) { return F(-1) / delta; });
  Dpphh->sum(F(1.0), "abij", "abij", [] (F delta) { return F(-1) / delta; });

  bool first(true);
  double accumulate_energy(0.0);
  auto compute_energy =
    [&Rph, &Rpphh, &Dph, &Dpphh, &Tph, &Tpphh, &Vpphh, &first]
    (const std::string s) -> double {

    auto ccdenergy( Tcc<TE>::template tensor<F>("ccdE") );
    auto ccxenergy( Tcc<TE>::template tensor<F>("ccxE") );
    auto ttdenergy( Tcc<TE>::template tensor<F>("ttdE") );
    auto ttxenergy( Tcc<TE>::template tensor<F>("ttxE") );
    COMPILE(
      (*ttdenergy)[""] <<= ( 2.0) * (*Rpphh)["abij"] * (*Tpphh)["abij"],
      (*ttxenergy)[""] <<= (-1.0) * (*Rpphh)["abij"] * (*Tpphh)["baij"],
      (*Rpphh)["abij"] <<= (*Rpphh)["abij"] * (*Dpphh)["abij"],
      (*Rph)["ai"]     <<= (*Rph)["ai"]     * (*Dph)["ai"],
      (*ccdenergy)[""] <<= ( 2.0) * (*Rpphh)["abij"] * (*Vpphh)["abij"],
      (*ccdenergy)[""]  += ( 2.0) * (*Rph)["ai"] * (*Rph)["bj"] * (*Vpphh)["abij"],
      (*ccxenergy)[""] <<= (-1.0) * (*Rpphh)["abij"] * (*Vpphh)["baij"],
      (*ccxenergy)[""]  += (-1.0) * (*Rph)["ai"] * (*Rph)["bj"] * (*Vpphh)["baij"]
    )->execute();
    F _ccd(ccdenergy->read()), _ccx(ccxenergy->read()),
      _ttd(ttdenergy->read()), _ttx(ttxenergy->read());
    if (first) {
      OUT() << _FORMAT("%15s       Method            CCTotal   CCDirect   CCExchange   T2Total    T2Direct   T2Exchange\n", " ");
      first = false;
    }

    OUT() << _FORMAT("Channel Energy %20s : %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n",
                      s.c_str(), real(_ccd+_ccx), real(_ccd), real(_ccx), real(_ttd+_ttx), real(_ttd), real(_ttx));
    return real(_ccd+_ccx);
  };



  auto r3mode( arguments->template getValue<std::string>("r3-mode", "ct-full"));
  auto t3t3mode( arguments->template getValue<std::string>("t3t3-mode", "full"));


  auto Vpppp(coulombIntegrals->get("pppp"));
  auto Vphph(coulombIntegrals->get("phph"));
  auto Vhhhh(coulombIntegrals->get("hhhh"));
  auto Vhhhp(coulombIntegrals->get("hhhp"));
  auto Vhphh(coulombIntegrals->get("hphh"));
  auto Vppph(coulombIntegrals->get("ppph"));
  auto Vhhpp(coulombIntegrals->get("hhpp"));
  auto Vpphp(coulombIntegrals->get("pphp"));
  auto Vphhh(coulombIntegrals->get("phhh"));
  auto Vhphp(coulombIntegrals->get("hphp"));
  auto Vhpph(coulombIntegrals->get("hpph"));
  auto Vphhp(coulombIntegrals->get("phhp"));
  auto Vphpp(coulombIntegrals->get("phpp"));
  auto Vhhph(coulombIntegrals->get("hhph"));
  auto Vhppp(coulombIntegrals->get("hppp"));

  auto Xabcd( Tcc<TE>::template tensor<F>("Xabcd") );

  COMPILE(
    (*Rph)["ai"] <<= 0.0 * (*Vpphh)["aaii"],
    (*Rpphh)["abij"] <<= (*Vpphh)["abij"]
  )->execute();
  accumulate_energy += compute_energy("MP2");

  COMPILE(
    (*Rph)["ai"] <<= 0.0 * (*Vpphh)["aaii"],
    (*Rpphh)["abij"] <<= (*Vhhhh)["klij"] * (*Tpphh)["abkl"]
  )->execute();
  accumulate_energy += compute_energy("hhl");

  COMPILE(
    (*Rph)["ai"] <<= 0.0 * (*Vpphh)["aaii"],
    (*Xabcd)["abcd"] <<= (1.0) * (*Vpppp)["abcd"],
    (*Xabcd)["abcd"]  += (-1.0) * (*Vphpp)["akcd"] * (*Tph)["bk"],
    (*Xabcd)["abcd"]  += (-1.0) * (*Vhppp)["kbcd"] * (*Tph)["ak"],
    (*Rpphh)["abij"] <<= (*Xabcd)["abcd"] * (*Tpphh)["cdij"],
    (*Rpphh)["abij"]  += (*Xabcd)["abcd"] * (*Tph)["ci"] * (*Tph)["dj"]
  )->execute();
  compute_energy("PPL (with singles)");

  COMPILE(
    (*Rph)["ai"] <<= 0.0 * (*Vpphh)["aaii"],
    (*Xabcd)["abcd"] <<= (1.0) * (*Vpppp)["abcd"],
    (*Rpphh)["abij"] <<= (*Xabcd)["abcd"] * (*Tpphh)["cdij"]
  )->execute();
  accumulate_energy += compute_energy("PPL (no singles)");

  COMPILE(
    (*Rph)["ai"] <<= 0.0 * (*Vpphh)["aaii"],
    (*Rpphh)["abij"] <<= (-1.0) * (*Vphph)["akci"] * (*Tpphh)["cbkj"],
    (*Rpphh)["abij"] += (*Rpphh)["baji"]
  )->execute();
  accumulate_energy += compute_energy("PHL (inside bubble)");

  COMPILE(
    (*Rph)["ai"] <<= 0.0 * (*Vpphh)["aaii"],
    (*Rpphh)["abij"] <<= (-1.0) * (*Vphph)["bkci"] * (*Tpphh)["ackj"],
    (*Rpphh)["abij"] += (*Rpphh)["baji"]
  )->execute();
  accumulate_energy += compute_energy("PHL (left-right)");

  COMPILE(
    (*Rph)["ai"] <<= 0.0 * (*Vpphh)["aaii"],
    (*Rpphh)["abij"] <<= (2.0) * (*Vphhp)["akic"] * (*Tpphh)["cbkj"],
    (*Rpphh)["abij"] += (*Rpphh)["baji"]
  )->execute();
  accumulate_energy += compute_energy("Ring");

  COMPILE(
    (*Rph)["ai"] <<= 0.0 * (*Vpphh)["aaii"],
    (*Rpphh)["abij"] <<= (-1.0) * (*Vphhp)["akic"] * (*Tpphh)["bckj"],
    (*Rpphh)["abij"] += (*Rpphh)["baji"]
  )->execute();
  accumulate_energy += compute_energy("xRing");

  COMPILE(
    (*Rph)["ai"] <<= 0.0 * (*Vpphh)["aaii"],
    (*Rpphh)["abij"] <<= (2.0) * (*Vhhpp)["ijab"] * (*Tpphh)["acik"] * (*Tpphh)["dblj"]
  )->execute();
  accumulate_energy += compute_energy("qu. Ring");


  OUT() << "Accumulated energy " << accumulate_energy << std::endl;

  auto result(New<MapNode>(SOURCE_LOCATION));
  result->setValue("allGood", 1);
  return result;


}

