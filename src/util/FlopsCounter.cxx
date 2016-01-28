#include <util/FlopsCounter.hpp>

using namespace cc4s;

FlopsCounter::FlopsCounter(
  int64_t *flops_, MPI_Comm comm_
): flops(flops_), comm(comm_) {
}

FlopsCounter::~FlopsCounter() {
  *flops = counter.count(comm);
}

