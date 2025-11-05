#include <StepCheck.hpp>

namespace cc4s {

static void check_perturbative_triples(Ptr<MapNode> &steps,
                                       Ptr<MapNode> &storage) {

  for (Natural<> j(0); j < steps->getSize(); ++j) {
    const Ptr<MapNode> step(steps->getMap(j));
    const std::string name = step->getSymbol("name");
    if (name != "PerturbativeTriples")
      continue;
    auto in = step->getMap("in");
    const bool step_is_last = steps->getSize() - 1 == j;
    // if step is last, it will delete the tensors always
    in->setValue<int>("deleteTensors", step_is_last ? 1 : 0);
  }
  return;
}

void StepCheck::check_all(Ptr<MapNode> &steps, Ptr<MapNode> &storage) {
  check_perturbative_triples(steps, storage);
}

} // namespace cc4s
