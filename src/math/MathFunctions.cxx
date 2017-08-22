#include <math/MathFunctions.hpp>

#include <util/Log.hpp>

int cc4s::permutationSign(
  const std::string &original, const std::string &permuted
) {
  // TODO: should this not report 0 if indices occur multiple times?
  int sign(+1);
  const int N(original.length());
  char prev;
  char next;
  std::string cycle(""); // cycle is not needed apart from length() and [0]
  std::string visited("");
  bool endOfCycle(false);
  // At most we will find length cycles
  for (int c(0); c < N; ++c) {
    cycle = "";
    endOfCycle = false;
    prev = original[c];
    next = permuted[c];
    LOG(1,"permutationSign") << "prev " << prev << " next " << next<< std::endl;
    if (visited.find(prev) != std::string::npos) {
      // If prev is already visited, skip it
      continue;
    }
    if (prev == next) {
      // Do not bother with fix points
      visited += prev;
      continue;
    }
    for (int i(0); i < N; ++i) {
      for (int j(0); j < N; ++j) {
        if (original[j] == next) {
          cycle += prev;
          sign = -sign;
          next = permuted[j];
          prev = original[j];
          if (next == cycle[0]) {
            cycle += prev;
            endOfCycle = true;
            continue;
          }
        }
      }
      if (endOfCycle) {
        LOG(1,"permutationSign") << "\tCycle = " << cycle << std::endl;
        break;
      }
    }
    visited += cycle;
    LOG(1,"permutationSign") << "visited  " << visited << std::endl;
  }
  LOG(1,"permutationSign") << "P ( " << original << " -> " << permuted << " ) "
            << "  =>  " << sign << std::endl;
  return sign;
}
