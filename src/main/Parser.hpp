#ifndef PARSER_DEFINED
#define PARSER_DEFINED

#include <Data.hpp>
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
    // TODO: parse from stream
    // TODO: maybe store source location info (file,line,column) in node
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
      switch (yamlNode.Type()) {
      case YAML::NodeType::Map:
        return parseMap(yamlNode);
      case YAML::NodeType::Sequence:
        return parseSequence(yamlNode);
      case YAML::NodeType::Scalar:
        return parseScalar(yamlNode);
      default:
        Assert(false, "unknown node type: " + yamlNode.Type());
      }
    }

    Ptr<MapNode> parseMap(const YAML::Node &yamlNode) {
      auto node(New<MapNode>());
      for (auto iterator: yamlNode) {
        node->get(iterator.first.as<std::string>()) = parseNode(iterator.second);
      }
      return node;
    }

    Ptr<MapNode> parseSequence(const YAML::Node &yamlNode) {
      auto node(New<MapNode>());
      size_t index(0);
      for (auto iterator: yamlNode) {
        LOG(1,"Parser") << "parsing sequence #" << index << std::endl;
        node->get(index++) = parseNode(iterator);
      }
      return node;
    }

    Ptr<Node> parseScalar(const YAML::Node &yamlNode) {
      // use explicitly node type, if given
      auto tag(yamlNode.Tag());
      if (tag.find(":str") != std::string::npos) {
        // by default, strings are symbol names
        return parseSymbol(yamlNode);
      } else if (tag.find(":text") != std::string::npos) {
        // literal strings have to be given the !!text tag
        return parseAtom<std::string>(yamlNode);
      } else if (tag.find(":int") != std::string::npos) {
        // default integers
        return parseAtom<int64_t>(yamlNode);
      } else if (
        tag.find(":float") != std::string::npos  ||
        tag.find(":real") != std::string::npos
      ) {
        // default reals
        return parseAtom<Real<>>(yamlNode);
      } else if (tag.find(":real64") != std::string::npos) {
        // explicit size reals
        return parseAtom<Real<64>>(yamlNode);
/*
      } else if (tag.find(":real128") != std::string::npos) {
        // TODO: 128 bit real parsing
        // explicity size reals
        return parseAtom<Real<128>>(yamlNode);
*/
      } else if (tag.find(":complex") != std::string::npos) {
        // default  complex
        return parseAtom<Complex<>>(yamlNode);
      } else if (tag.find(":complex64") != std::string::npos) {
        // explicit size complex
        return parseAtom<Complex<64>>(yamlNode);
/*
      } else if (tag.find(":complex128") != std::string::npos) {
        // TODO: 128 bit real parsing
        // explicit size complex
        return parseAtom<Complex<128>>(yamlNode);
*/
      } else if (tag == "?" ){
        return parseImplicitTypeScalar(yamlNode);
      } else {
        auto value(yamlNode.Scalar());
        Assert(
          false, "unsupported type tag " + tag + " for value " + value
        );
      }
    }

    Ptr<Node> parseImplicitTypeScalar(const YAML::Node &yamlNode) {
      // implicitly deduce node type from value
      auto value(yamlNode.Scalar());
      if (value.size() == 0) {
        // empty string can only be text
        return parseAtom<std::string>(yamlNode);
      } else if (isdigit(value[0]) || value[0] == '-' || value[0] == '+') {
        return parseNumber(yamlNode);
      } else {
        return parseSymbol(yamlNode);
      }
    }

    Ptr<Node> parseNumber(const YAML::Node &yamlNode) {
      if (yamlNode.Scalar().find('(') != std::string::npos) {
        return parseAtom<Complex<>>(yamlNode);
      } else if (yamlNode.Scalar().find('.') != std::string::npos) {
        return parseAtom<Real<>>(yamlNode);
      } else {
        return parseAtom<int64_t>(yamlNode);
      }
    }

    Ptr<Node> parseSymbol(const YAML::Node &yamlNode) {
      auto value(yamlNode.Scalar());
      if (value.size() == 0) {
        // empty string can only be text
        return parseAtom<std::string>(yamlNode);
      } else if (isalpha(value[0])) {
        // first character is alphabetic
        size_t i(1);
        // rest must be alphabetic or digit
        while ((isalpha(value[i]) || isdigit(value[i])) && i<value.size()) ++i;
        if (i == value.size()) {
          // then it's a symbol
          return New<SymbolNode>(yamlNode.Scalar());
        }
      }
      // otherwise it's text
      return parseAtom<std::string>(yamlNode);
    }

    template <typename AtomicType>
    Ptr<AtomicNode<AtomicType>> parseAtom(const YAML::Node &yamlNode) {
      return New<AtomicNode<AtomicType>>(yamlNode.as<AtomicType>());
    }

    std::string fileName;
  };
}

// provide conversion to and from complex
// FIXME: doesn't work properly
namespace YAML {
template<>
struct convert<cc4s::Complex<64>> {
  static Node encode(const cc4s::Complex<64>& rhs) {
    std::stringstream stream;
    stream << rhs;
    Node node;
    node = stream.str();
    return node;
  }

  static bool decode(const Node& node, cc4s::Complex<64>& rhs) {
    std::stringstream stream(node.Scalar());
    stream >> rhs;
    return true;
  }
};
}
#endif

