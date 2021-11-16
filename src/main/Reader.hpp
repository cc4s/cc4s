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

#ifndef READER_DEFINED
#define READER_DEFINED

#include <Node.hpp>
#include <Parser.hpp>
#include <Cc4s.hpp>
#include <util/SharedPointer.hpp>
#include <string>

namespace cc4s {
  /**
   * \brief Reader for reading entire node states from files,
   * including tensor data.
   */
  class Reader {
  public:
    Reader(
      const std::string &pathFileName_
    ): pathFileName(pathFileName_) {
    }
    Ptr<Node> read() {
      auto dotPosition(pathFileName.rfind('.'));
      ASSERT(
        dotPosition != std::string::npos,
        "'fileName' must contain an extension e.g. '.yaml'"
      );
      // if pathFileName contains '/' change directory
      char currentDirectory[PATH_MAX];
      getcwd(currentDirectory, sizeof(currentDirectory));
      // TODO: support other filesystems
      std::string fileName, baseName;
      auto slashPosition(pathFileName.rfind('/'));
      if (slashPosition != std::string::npos) {
        auto fileDirectory(pathFileName.substr(0, slashPosition));
        chdir(fileDirectory.c_str());
        fileName = pathFileName.substr(slashPosition+1);
        baseName = pathFileName.substr(slashPosition+1, dotPosition-slashPosition-1);
      } else {
        fileName = pathFileName;
        baseName = pathFileName.substr(0, dotPosition);
      }

      // read persistent tree from given file
      auto persistentNode(Parser(fileName).parse());

      // translate from persistent node tree into working node tree.
      // e.g. tensors are read into AtomicNode<Ptr<Tensor<F,TE>> etc.
      // Contained tensor data are read as side effects.
      auto node(read(persistentNode, baseName));

      // restore to original directory
      chdir(currentDirectory);

      return node;
    }

    typedef std::function<
      Ptr<Node>(const Ptr<MapNode> &node, const std::string &nodePath)
    > ReadFunction;
    static int registerReadFunction(
      const std::string &name, ReadFunction readFunction
    ) {
      readFunctions[name] = readFunction;
      return 1;
    }

  protected:
    static std::map<std::string, ReadFunction> readFunctions;

    // translate into working node tree
    Ptr<Node> read(const Ptr<Node> &node, const std::string &nodePath) {
      auto mapNode(node->toPtr<MapNode>());
      if (mapNode) return readMap(mapNode, nodePath);
      // other nodes are already in working form
      return node;
    }

    Ptr<Node> readMap(
      const Ptr<MapNode> &mapNode,
      const std::string &nodePath
    ) {
      auto readFunction(readFunctions.end());
      // check if map has key "type" and whether there is a read function
      // registered for that type
      if (mapNode->get("type")) {
        readFunction = readFunctions.find(
          mapNode->getValue<std::string>("type")
        );
      }
      if (readFunction != readFunctions.end()) {
        // if yes, call respective handler read function to translate into
        // node of working tree
        return readFunction->second(mapNode, nodePath);
      } else {
        // otherwise: translate sub nodes
        auto workingMapNode(New<MapNode>(mapNode->sourceLocation));
        for (auto key: mapNode->getKeys()) {
          auto subNode(mapNode->get(key));
          // NOTE that there may be nullptr keys, if keys have been searched for
          if (subNode) {
            workingMapNode->get(key) = read(subNode, nodePath + "." + key);
          }
        }
        return workingMapNode;
      }
    }

    std::string pathFileName;
  };
}

#endif

