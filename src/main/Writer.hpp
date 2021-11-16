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

#ifndef WRITER_DEFINED
#define WRITER_DEFINED

#include <Node.hpp>
#include <Emitter.hpp>
#include <Cc4s.hpp>
#include <util/SharedPointer.hpp>
#include <string>

namespace cc4s {
  /**
   * \brief Writer for writing entire node states persistently to files,
   * including tensor data.
   * In contrast, the Emitter only writes memory references of tensors,
   * which can only be used during the execution of cc4s.
   */
  class Writer {
  public:
    Writer(
      const std::string &fileName_, const bool useBinary_ = false
    ): fileName(fileName_), useBinary(useBinary_) {
    }
    Ptr<Node> write(const Ptr<Node> &node) {
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
      Emitter(fileName).emit(writtenNode);

      return writtenNode;
    }

    typedef std::function<
      Ptr<Node>(
        const Ptr<Node> &node, const std::string &nodePath, const bool useBinary
      )
    > WriteFunction;
    static int registerWriteFunction(
      const std::string &name, WriteFunction writeFunction
    ) {
      writeFunctions[name] = writeFunction;
      return 1;
    }

  protected:
    static std::map<std::string, WriteFunction> writeFunctions;

    // translate into persistent node tree
    Ptr<Node> write(const Ptr<Node> &node, const std::string &nodePath) {
      auto mapNode(node->toPtr<MapNode>());
      if (mapNode) return writeMap(mapNode, nodePath);
      auto symbolNode(node->toPtr<SymbolNode>());
      if (symbolNode) {
        // symbols are already in persistent form
        return symbolNode;
      } else {
        return writeAtom(node, nodePath);
      }
    }

    Ptr<MapNode> writeMap(
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

    Ptr<Node> writeAtom(
      const Ptr<Node> &node, const std::string &nodePath
    ) {
      for (auto writeFunction: writeFunctions) {
        // try type-specific write function
        auto writtenNode(writeFunction.second(node, nodePath, useBinary));
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

    std::string fileName;
    bool useBinary;
  };
}

#endif

