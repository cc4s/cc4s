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
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <math/TensorUnion.hpp>
#include <math/MathFunctions.hpp>
#include <atrip.hpp>
#include <atrip/Debug.hpp>

using namespace cc4s;
ALGORITHM_REGISTRAR_DEFINITION(PerturbativeTriples)

#define Q(...) #__VA_ARGS__
#define QUOTE(...) Q(__VA_ARGS__)

Ptr<MapNode> PerturbativeTriples::run(const Ptr<MapNode> &arguments) {
  using TE = DefaultTensorEngine;

  auto result(New<MapNode>(SOURCE_LOCATION));

  atrip::Atrip::init();
  atrip::Atrip::Input in;

#define __V__(_idx)                                         \
    ([&arguments]() {                                       \
      COMPILE(arguments                                     \
        ->getMap("coulombIntegrals")                        \
        ->getMap("slices")                                  \
        ->getPtr<TensorExpression<Real<>,TE>>(#_idx)        \
        )->execute();                                       \
      return                                                \
          &(arguments                                       \
            ->getMap("coulombIntegrals")                    \
            ->getMap("slices")                              \
            ->getPtr<TensorExpression<Real<>,TE>>(#_idx)    \
            ->evaluate()                                    \
            ->getMachineTensor()                            \
            ->tensor);                                      \
    })()
#define __T__(_idx) \
   &(arguments                                                    \
      ->getPtr<TensorUnion<Real<>,TE>>("amplitudes")              \
      ->get(_idx)                                                 \
      ->evaluate()                                                \
      ->getMachineTensor()                                        \
      ->tensor)
#define __eps__(_idx)                                      \
   &(arguments                                             \
      ->getMap("slicedEigenEnergies")                      \
      ->getMap("slices")                                   \
      ->getPtr<TensorExpression<Real<>,TE>>(#_idx)         \
      ->evaluate()                                         \
      ->getMachineTensor()                                 \
      ->tensor)

  // this is a hack so that a CTF::World gets created for sure
  CTF::World _w(MPI_COMM_WORLD);
  in
    // setup tensors
    .with_epsilon_i(__eps__(h))
    .with_epsilon_a(__eps__(p))
    .with_Tai(__T__(0))
    .with_Tabij(__T__(1))
    .with_Vabij(__V__(pphh))
    .with_Vijka(__V__(hhhp))
    .with_Vabci(__V__(ppph))
    // some options
    .with_barrier(false)
    .with_percentageMod(10)
    ;

  const double No = double((__V__(pphh))->lens[2])
             , Nv = double((__V__(pphh))->lens[0])
             , doublesFlops = No * No * No
                            * (No + Nv)
                            * 2.0
                            * 6.0
                            / 1e9
                            ;
  //OUT() << "doublesFlops " << doublesFlops << std::endl;
  auto const rank = Cc4s::world->getRank();
  bool firstHeaderPrinted = false;
  atrip::registerIterationDescriptor
    ([doublesFlops, &firstHeaderPrinted, rank]
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
              d.currentElapsedTime,
              doublesFlops / d.currentElapsedTime);
      OUT() << out << "\n";
    });


#undef __V__
#undef __T__
#undef __eps__

  auto out = atrip::Atrip::run(in);
  OUT() << "(T) correlation energy: "
        << std::setprecision(15) << std::setw(23)
        << out.energy << std::endl;
  OUT() << "(T*) correlation energy: " << "TODO" << std::endl;

  auto energy(New<MapNode>(SOURCE_LOCATION));
  energy->setValue("triples", real(out.energy));
  result->get("energy") = energy;

  return result;
}
