#ifndef ___STEPCHECK_HPP
#define ___STEPCHECK_HPP

#include <Node.hpp>

namespace cc4s {

class StepCheck {
public:
  static bool PerturbativeTriplesIsLast;

  static void check_all(Ptr<MapNode> &steps, Ptr<MapNode> &storage);
};

} // namespace cc4s

#endif
