/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include "Emitter.hpp"

#include <fstream>
#include <string>

using namespace cc4s;

int Emitter::rank(-1);
PTR(std::ofstream) Emitter::yamlFile;
PTR(YAML::Emitter) Emitter::yamlEmitter;

void Emitter::setRank(int const rank_) {
  rank = rank_;
}

int Emitter::getRank() {
  return rank;
}

YAML::Emitter &Emitter::getEmitter() {
  if (!yamlFile) {
    yamlFile = NEW(
      std::ofstream, "cc4s.yaml", std::ofstream::out | std::ofstream::trunc
    );
  }
  if (!yamlEmitter) {
    yamlEmitter = NEW(YAML::Emitter, *yamlFile);
  }
  return *yamlEmitter;
}

