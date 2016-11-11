#include <util/BlacsWorld.hpp>

#include <extern/Blacs.hpp>
#include <extern/ScaLapack.hpp>
#include <util/Log.hpp>

using namespace cc4s;

BlacsWorld::BlacsWorld(int rank_, int processes, int processRows): rank(rank_) {
  lens[0] = processRows > 0 ?
    processRows : static_cast<int>(std::sqrt(processes));
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

