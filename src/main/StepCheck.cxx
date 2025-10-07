#include <StepCheck.hpp>

namespace cc4s {

bool StepCheck::PerturbativeTriplesIsLast = true;
static void check_perturbative_triples(Ptr<MapNode> &steps,
                                       Ptr<MapNode> &storage) {

  for (Natural<> j(0); j < steps->getSize(); ++j) {
    const Ptr<MapNode> step(steps->getMap(j));
    const std::string name = step->getSymbol("name");
    if (name != "PerturbativeTriples")
      continue;
    StepCheck::PerturbativeTriplesIsLast = steps->getSize() - 1 == j;
    return;
  }
}

void StepCheck::check_all(Ptr<MapNode> &steps, Ptr<MapNode> &storage) {
  check_perturbative_triples(steps, storage);
}

} // namespace cc4s
