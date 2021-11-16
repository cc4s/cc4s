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

#ifndef EMITTER_DEFINED
#define EMITTER_DEFINED

#include <Node.hpp>
#include <Cc4s.hpp>
#include <util/SharedPointer.hpp>
#include <string>
#include <yaml-cpp/yaml.h>

namespace cc4s {
  /**
   * \brief Emitter for emitting the data nodes to a yaml file
   */
  class Emitter {
  public:
    Emitter(const std::string &fileName): yamlEmitter(stream) {
      if (Cc4s::world->getRank() == 0) stream.open(fileName);
    }
    void emit(const Ptr<Node> &node) {
      auto mapNode(node->toPtr<MapNode>());
      if (mapNode) return emitMap(mapNode);
      auto symbolNode(node->toPtr<SymbolNode>());
      if (symbolNode) {
        emitSymbol(symbolNode);
      } else {
        emitAtom(node);
      }
    }
  protected:
    void emitSymbol(Ptr<SymbolNode> symbolNode) {
      yamlEmitter << symbolNode->value;
    }

    void emitAtom(Ptr<Node> node) {
      // TODO: emit comments in node structure
      // TODO: emit with explicit type tags
      yamlEmitter << node->toString();
    }

    void emitMap(Ptr<MapNode> mapNode) {
      yamlEmitter << YAML::BeginMap;
      for (auto key: mapNode->getKeys()) {
        auto valueNode(mapNode->get(key));
        // NOTE that there may be nullptr keys, if keys have been searched for
        if (valueNode) {
          yamlEmitter << YAML::Key << key;
          yamlEmitter << YAML::Value;
          emit(valueNode);
        }
      }
      yamlEmitter << YAML::EndMap;
    }

    std::ofstream stream;
    YAML::Emitter yamlEmitter;
  };
}

#endif

