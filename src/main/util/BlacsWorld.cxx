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

#include <util/BlacsWorld.hpp>

#include <extern/Blacs.hpp>
#include <extern/ScaLapack.hpp>
#include <util/Log.hpp>

using namespace cc4s;

// only if we have cblacs and scalapack in the project
#if defined(HAVE_CBLACS) && defined(HAVE_SCALAPACK)

BlacsWorld::BlacsWorld(int rank_, int processes, int processRows): rank(rank_) {
  lens[0] = processRows > 0 ?
    processRows : static_cast<int>(sqrt(processes));
  lens[1] = processes / lens[0];
  while (lens[0] * lens[1] != processes) {
    --lens[0];
    lens[1] = processes / lens[0];
  }
  Cblacs_get(-1, 0, &context);
  Cblacs_gridinit(&context, "ColumnMajor", lens[0], lens[1]);
  Cblacs_gridinfo(
    context, &lens[0], &lens[1], &firstElement[0], &firstElement[1]
  );
}

BlacsWorld::~BlacsWorld() {
  Cblacs_gridexit(context);
}

void BlacsWorld::barrier() {
  Cblacs_barrier(context, "All");
}

#endif
