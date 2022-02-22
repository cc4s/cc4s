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

#include <algorithms/PerturbativeTriples.hpp>
#include <algorithms/PerturbativeTriplesStar.hpp>
#include <algorithms/PerturbativeTriplesReference.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <TensorSet.hpp>
#include <MathFunctions.hpp>
#include <iomanip>
#include <atrip.hpp>
#include <atrip/Debug.hpp>

using namespace cc4s;
ALGORITHM_REGISTRAR_DEFINITION(PerturbativeTriples)

template <typename F, typename TE>
Ptr<MapNode> runAtrip(const Ptr<MapNode> &arguments) {
  {
    auto amplitudes = arguments->getPtr<TensorSet<F,TE>>("amplitudes");
    if (!amplitudes) return nullptr;
  }

  auto result(New<MapNode>(SOURCE_LOCATION));

  const auto
    distribution
      = arguments->getValue<std::string>("tuplesDistribution", "group") == "group"
      ? atrip::Atrip::Input<F>::TuplesDistribution::GROUP_AND_SORT
      : atrip::Atrip::Input<F>::TuplesDistribution::NAIVE
      ;

  const bool
    rankRoundRobin
      = arguments->getValue<std::string>("tuplesRoundRobin", "node") == "node"
      ? false
      : true
      ;

  atrip::Atrip::init();
  atrip::Atrip::Input<F> in;

#define __V__(_idx)                                            \
    ([&arguments]() {                                          \
      COMPILE(arguments                                        \
        ->getPtr<TensorSet<F,TE>>("coulombIntegrals")          \
        ->get(#_idx)                                           \
        )->execute();                                          \
      return                                                   \
          &(arguments                                          \
            ->getPtr<TensorSet<F,TE>>("coulombIntegrals")      \
            ->get(#_idx)                                       \
            ->evaluate()                                       \
            ->getMachineTensor()                               \
            ->tensor);                                         \
    })()
#define __T__(_idx)                                          \
   &(arguments                                               \
      ->getPtr<TensorSet<F,TE>>("amplitudes")                \
      ->get(_idx)                                            \
      ->evaluate()                                           \
      ->getMachineTensor()                                   \
      ->tensor)
#define __eps__(_idx)                                        \
   &(arguments                                               \
      ->getPtr<TensorSet<Real<>,TE>>("slicedEigenEnergies")  \
      ->get(#_idx)                                           \
      ->evaluate()                                           \
      ->getMachineTensor()                                   \
      ->tensor)


  CTF::Tensor<F>
      epsi(1, (__eps__(h))->lens, (__eps__(h))->sym)
    , epsa(1, (__eps__(p))->lens, (__eps__(p))->sym)
    ;
  const auto toComplex
    = CTF::Transform<double, F>([](double d, F &f) { f = d; });
  toComplex((*__eps__(h))["i"], epsi["i"]);
  toComplex((*__eps__(p))["a"], epsa["a"]);

  // this is a hack so that a CTF::World gets created for sure
  CTF::World _w(MPI_COMM_WORLD);
  in
    // setup tensors
    .with_epsilon_i(&epsi)
    .with_epsilon_a(&epsa)
    .with_Tai(__T__("ph"))
    .with_Tabij(__T__("pphh"))
    .with_Vabij(__V__(pphh))
    .with_Vijka(__V__(hhhp))
    .with_Vabci(__V__(ppph))
    // some options
    .with_barrier(false)
    .with_percentageMod(10)
    .with_tuplesDistribution(distribution)
    .with_rankRoundRobin(rankRoundRobin)
    ;

  const double _No = *(__eps__(h))->lens
             , _Nv = *(__eps__(p))->lens
             , doublesFlops = _No * _No * _No
                            * (_No + _Nv)
                            * double(sizeof(F) / sizeof(double))
                            * 2.0
                            * 6.0
                            / 1.0e9
                            ;


  double lastElapsedTime = 0;
  auto const rank = Cc4s::world->getRank();
  bool firstHeaderPrinted = false;
  atrip::registerIterationDescriptor
    ([doublesFlops, &firstHeaderPrinted, rank, &lastElapsedTime]
       (atrip::IterationDescription const& d) {
      const char
        *fmt_header = "%-13s%-10s%-13s",
        *fmt_nums = "%-13.0f%-10.0f%-13.3f";
      char out[256];
      if (!firstHeaderPrinted) {
        sprintf(out, fmt_header, "Progress(%)", "time(s)", "GFLOP/s");
        firstHeaderPrinted = true;
        OUT() << out << "\n";
      }
      sprintf(out, fmt_nums,
              double(d.currentIteration) / double(d.totalIterations) * 100,
              (d.currentElapsedTime - lastElapsedTime),
              d.currentIteration * doublesFlops / d.currentElapsedTime);
      lastElapsedTime = d.currentElapsedTime;
      OUT() << out << "\n";
    });

#undef __V__
#undef __T__
#undef __eps__

  auto out = atrip::Atrip::run<F>(in);
  OUT() << "(T) correlation energy: "
        << std::setprecision(15) << std::setw(23)
        << out.energy << std::endl;

  /* DUPLICATION WARNING:
   * --------------------
   *
   * If you change something in the output of this algorithm
   * be sure to change accordingly the output of the
   * PerturbativeTriplesReference, as they should be
   * symmetric.
   *
   */

  auto energy(New<MapNode>(SOURCE_LOCATION));
  energy->setValue("correlation", std::real(out.energy));

  if (arguments->isGiven("mp2PairEnergies")) {
    const F
      triples_star = computeCssdPtStar<F, TE>(arguments, out.energy);
    OUT() << "(T*)-Bsie energy correction: "
          << std::setprecision(15) << std::setw(23)
          << std::real(triples_star) << std::endl;
    energy->setValue("starCorrection", std::real(triples_star));
  }

  using TSr = TensorSet<Real<>, TE>;
  energy
    ->setValue("unit",
               arguments
                ->getPtr<TSr>("slicedEigenEnergies")
                ->get("h")
                ->inspect()
                ->getUnit());

  result->get("energy") = energy;

  return result;
}

template <typename F, typename TE>
Ptr<MapNode> atripDryRun(const Ptr<MapNode> &arguments) {
  auto result(New<MapNode>(SOURCE_LOCATION));
  auto eps(arguments->getPtr<TensorSet<Real<>,TE>>("slicedEigenEnergies"));
  auto epsh(eps->get("h"));
  auto epsp(eps->get("p"));
  const
  size_t
      no = epsh->inspect()->getLen(0)
    , nv = epsp->inspect()->getLen(0)
    , nranks = Cc4s::getProcessesCount()
    , f = sizeof(F)
    , n_tuples = nv * (nv + 1) * (nv + 2) / 6 - nv
    , atrip_memory
        = /* tuples_memory */ 3 * sizeof(size_t) * n_tuples
          //
          // one dimensional slices (all ranks)
          //
        + /* taphh */ f * nranks * 6 * nv * no * no
        + /* hhha  */ f * nranks * 6 * no * no * no
          //
          // two dimensional slices (all ranks)
          //
        + /* abph  */ f * nranks * 12 * nv * no
        + /* abhh  */ f * nranks *  6 * no * no
        + /* tabhh */ f * nranks *  6 * no * no
          //
          // distributed sources (all ranks)
          //
        + /* tpphh */ f * nv * nv * no * no
        + /* vhhhp */ f * no * no * no * nv
        + /* vppph */ f * nv * nv * nv * no
        + /* vpphh */ f * nv * nv * no * no
        + /* tpphh2 */ f * nv * nv * no * no
          //
          // tensors in every rank
          //
        + /* tijk */ f * nranks * no * no * no
        + /* zijk */ f * nranks * no * no * no
        + /* epsp */ f * nranks * (no + nv)
        + /* tai  */ f * nranks * no * nv
    ;
  DryMemory::allocate(atrip_memory, SOURCE_LOCATION);
  DryMemory::free(atrip_memory);
  return result;
}

Ptr<MapNode> PerturbativeTriples::run(const Ptr<MapNode> &arguments) {

  Ptr<MapNode> result;

  if (Cc4s::dryRun) {
    using TE = DefaultDryTensorEngine;
    if (arguments->getPtr<TensorSet<Real<>,TE>>("amplitudes") != nullptr) {
      result = atripDryRun<Real<>, TE>(arguments);
    } else {
      result = atripDryRun<Complex<>, TE>(arguments);
    }
  } else {
    using TE = DefaultTensorEngine;
    (
      result = runAtrip<Real<>,TE>(arguments)
    ) || (
      result = runAtrip<Complex<>,TE>(arguments)
    );
  }

  return result;
}

