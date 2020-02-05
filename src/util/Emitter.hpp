/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef EMITTER_DEFINED
#define EMITTER_DEFINED

#include <util/SharedPointer.hpp>

#include <string>
#include <yaml-cpp/yaml.h>

namespace cc4s {
  /**
   * \brief Class with static members offering control over yaml emitting.
   * Entries are emitted with the macro EMIT.
   */
  class Emitter {
  public:
    static void setRank(const int rank);
    static int getRank();
    static YAML::Emitter &getEmitter();

  protected:
    static int rank;
    static PTR(std::ofstream) yamlFile;
    static PTR(YAML::Emitter) yamlEmitter;
  };
}

#define EMIT(...) \
  if (cc4s::Emitter::getRank() != 0) { \
  } else cc4s::Emitter::getEmitter()

#endif

