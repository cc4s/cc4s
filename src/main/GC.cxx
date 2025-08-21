#include <map>
#include <string>
#include <vector>

#include <GC.hpp>
#include <Log.hpp>

namespace cc4s {

void collect_garbage(size_t const i, Ptr<MapNode> &steps,
                     Ptr<MapNode> &storage) {

  std::map<std::string, bool> should_collect;
  // take a look at all the storage and see if there is
  // some algorithms coming that need any of the data
  // in the storage
  // Here we only work with the data name
  for (const auto &data_name : storage->getKeys()) {
    should_collect[data_name] = true;
    for (Natural<> j(i); j < steps->getSize(); ++j) {
      const Ptr<MapNode> step(steps->getMap(j));
      const Ptr<MapNode> inputArguments(step->getMap("in"));
      for (const auto &in_key_name : inputArguments->getKeys()) {
        const Ptr<Node> value_node = inputArguments->get(in_key_name);
        const std::string string_node = value_node->toString();
        if (string_node == data_name) {
          should_collect[data_name] = false;
        }
      }
    }
  }

  for (auto const &_should_collect_pair : should_collect) {
    const std::string data_name = _should_collect_pair.first;
    const bool should_collect_p = _should_collect_pair.second;
    if (!should_collect_p)
      continue;
    LOG() << "GC: garbage collecting " << data_name << "  (#"
          << storage->get(data_name).use_count() << ")"
          << "  :" << storage->get(data_name).get() << std::endl;

    for (Natural<> j(0); j < steps->getSize(); ++j) {
      Ptr<MapNode> step(steps->getMap(j));
      Ptr<MapNode> inputArguments(step->getMap("in"));
      Ptr<MapNode> outputArguments = step->isGiven("out")
                                         ? step->getMap("out")
                                         : New<MapNode>(step->sourceLocation);

      // remove in all input arguments
      for (auto args : {inputArguments, outputArguments})
        for (const auto &key_name : args->getKeys()) {
          Ptr<Node> *value_node = &args->get(key_name);
          if (value_node->get() == storage->get(data_name).get()) {
            value_node->reset();
            LOG() << "GC: Removed " << data_name << " in step " << j << "  (#"
                  << storage->get(data_name).use_count() << ")"
                  << "  :" << storage->get(data_name).get() << std::endl;
          }
        }
    }

    // delete in storage
    for (const auto &key_name : storage->getKeys()) {
      if (key_name == data_name) {
        storage->get(key_name).reset();
        LOG() << "GC: Removed " << data_name << " ("
              << storage->get(key_name).use_count() << ") in storage"
              << "  (#" << storage->get(data_name).use_count() << ")"
              << "  :" << storage->get(data_name).get() << std::endl;
      }
    }
  }
}

} // namespace cc4s
