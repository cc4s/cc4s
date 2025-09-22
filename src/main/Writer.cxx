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

#include <Writer.hpp>

namespace cc4s {

  std::map<std::string, Writer::WriteFunction> Writer::writeFunctions;

  int Writer::registerWriteFunction(const std::string &name,
                                    Writer::WriteFunction writeFunction) {
    writeFunctions[name] = writeFunction;
    return 1;
  }

  Ptr<Node> Writer::write(const Ptr<Node> &node) {
    // if fileName contains '/' change directory
    char currentDirectory[PATH_MAX];
    getcwd(currentDirectory, sizeof(currentDirectory));
    // TODO: support other filesystems
    auto slashPosition(fileName.rfind('/'));
    if (slashPosition != std::string::npos) {
      auto fileDirectory(fileName.substr(0, slashPosition));
      chdir(fileDirectory.c_str());
    }
    auto dotPosition(fileName.rfind('.'));
    ASSERT(
           dotPosition != std::string::npos,
           "'fileName' must contain an extension e.g. '.yaml'"
           );
    // base name needed for tensor data
    auto baseName(
                  fileName.substr(slashPosition+1, dotPosition-slashPosition-1)
                  );

    // translate into persistent node tree. Contained tensor data are
    // written as side effects.
    auto writtenNode(write(node, baseName));

    // restore to original directory
    chdir(currentDirectory);

    // finally, write persistent tree to given file
    if (Cc4s::dryRun)  fileName += "_";
    Emitter(fileName).emit(writtenNode);
    if (Cc4s::dryRun) remove(fileName.c_str());

    return writtenNode;
  }


  // translate into persistent node tree
  Ptr<Node> Writer::write(const Ptr<Node> &node, const std::string &nodePath) {
    auto mapNode(node->toPtr<MapNode>());
    if (mapNode) {
      return writeMap(mapNode, nodePath);
    }
    auto symbolNode(node->toPtr<SymbolNode>());
    if (symbolNode) {
      // symbols are already in persistent form
      ASSERT_LOCATION(0, std::string("unassigned object ") + symbolNode->value,
                      symbolNode->sourceLocation);
      return symbolNode;
    } else {
      return writeAtom(node, nodePath);
    }
  }

  Ptr<MapNode> Writer::writeMap(
                                const Ptr<MapNode> &mapNode,
                                const std::string &nodePath
                                ) {
    auto writtenMapNode(New<MapNode>(mapNode->sourceLocation));
    for (auto key: mapNode->getKeys()) {
      auto subNode(mapNode->get(key));
      // NOTE that there may be nullptr keys, if keys have been searched for
      if (subNode) {
        writtenMapNode->get(key) = write(subNode, nodePath + "." + key);
      }
    }
    return writtenMapNode;
  }


  Ptr<Node> Writer::writeAtom(
                              const Ptr<Node> &node, const std::string &nodePath
                              ) {
    for (auto writeFunction: writeFunctions) {
      // try type-specific write function
      auto writtenNode(writeFunction.second(node, nodePath, flagsIO));
      if (writtenNode) {
        // if successful, enter respective type
        auto writtenNodeMap(writtenNode->toPtr<MapNode>());
        writtenNodeMap->setSymbol("type", writeFunction.first);
        return writtenNode;
      }
    }
    // other types of atomic nodes can be serialized with to_string
    return node;
  }

}  // namespace cc4s
