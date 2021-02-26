/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef READER_DEFINED
#define READER_DEFINED

#include <Data.hpp>
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
      auto mapNode(node->toMap());
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

