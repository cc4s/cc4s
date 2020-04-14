#ifndef PARSER_DEFINED
#define PARSER_DEFINED

#include <Data.hpp>
#include <util/LineNumberStream.hpp>
#include <util/SharedPointer.hpp>
#include <yaml-cpp/yaml.h>
#include <string>
#include <istream>

namespace cc4s {
  /**
   * \brief Parser for cc4s files specifying the calculation plan, i.e.
   * which algorithms to use in which order.
   */
  class Parser {
  public:
    /**
     * \brief Creates a new interpreter for a cc4s yaml file of the given name.
     * Upon creation the file will be openend but not yet read.
     */
/*
    Parser(const std::string &fileName): stream(New<std::ifstream>(fileName)) {
      if (!std::dynamic_pointer_cast<std::ifstream>(stream)->is_open()) {
        Assert(false, "Failed to open file " + fileName);
      }
    }
*/
    Parser(const std::string &fileName_): fileName(fileName_) {
    }

    /**
     * \brief Parses the cc4s yaml data contained in the stream.
     * This method must be called with the same stream content on all processes.
     */
    Ptr<Node> parse() {
      LOG(1,"Parser") << "Parsing file " << fileName << std::endl;
      YAML::Node yamlNode(YAML::LoadFile(fileName));
      return parseNode(yamlNode);
    }

  protected:
    Ptr<Node> parseNode(const YAML::Node &yamlNode) {
      LOG(1,"Parser") << "node type: " << yamlNode.Type() << " " << YAML::NodeType::Scalar << std::endl;
      LOG(1,"Parser") << "node tag: " << yamlNode.Tag() << std::endl;
      switch (yamlNode.Type()) {
      case YAML::NodeType::Map:
        return parseMap(yamlNode);
      case YAML::NodeType::Sequence:
        return parseSequence(yamlNode);
      case YAML::NodeType::Scalar:
        if (yamlNode.Tag().find("str") != std::string::npos) {
          // by default, strings are symbol names
          return parseSymbol(yamlNode);
        } else if (yamlNode.Tag().find("text") != std::string::npos) {
          // literal strings have to be given the !!text tag
          return parseAtom<std::string>(yamlNode);
        } else if (yamlNode.Tag().find("int") != std::string::npos) {
          // default integers
          return parseAtom<int64_t>(yamlNode);
        } else if (yamlNode.Tag().find("float") != std::string::npos) {
          // default reals
          return parseAtom<Real<>>(yamlNode);
        } else if (yamlNode.Tag().find("real64") != std::string::npos) {
          // explicity size reals
          return parseAtom<Real<64>>(yamlNode);
        } else if (yamlNode.Tag().find("real128") != std::string::npos) {
          // explicity size reals
          // TODO: 128 bit real parsing
          return nullptr;
//          return parseAtom<Real<128>>(yamlNode);
        } else {
          return parseSymbol(yamlNode);
        }
      default:
        Assert(false, "unknown node type: " + yamlNode.Type());
      }
    }

    Ptr<MapNode> parseMap(const YAML::Node &yamlNode) {
      auto node(New<MapNode>());
      for (auto iterator: yamlNode) {
        LOG(1,"Parser") << "map key " << iterator.first.as<std::string>() << std::endl;
        node->get(iterator.first.as<std::string>()) = parseNode(iterator.second);
      }
      return node;
    }

    Ptr<MapNode> parseSequence(const YAML::Node &yamlNode) {
      auto node(New<MapNode>());
      size_t index(0);
      for (auto iterator: yamlNode) {
        LOG(1,"Parser") << "sequence number " << index << std::endl;
        node->get(index++) = parseNode(iterator);
      }
      return node;
    }

    Ptr<SymbolNode> parseSymbol(const YAML::Node &yamlNode) {
      // TODO: check validity of symbol name
      // TODO: create text node if not a valid symbol name
      return New<SymbolNode>(yamlNode.Scalar());
    }

    template <typename AtomicType>
    Ptr<AtomicNode<AtomicType>> parseAtom(const YAML::Node &yamlNode) {
      return New<AtomicNode<AtomicType>>(yamlNode.as<AtomicType>());
    }

    Ptr<std::istream> stream;
    std::string fileName;
  };
}

#endif

