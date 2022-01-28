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

#ifndef PARSER_DEFINED
#define PARSER_DEFINED

#include <Node.hpp>
#include <Log.hpp>
#include <SharedPointer.hpp>

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
      LOG() << "Parsing file " << fileName << std::endl;
      try {
        YAML::Node yamlNode(YAML::LoadFile(fileName));
        return parseNode(yamlNode);
      } catch (const YAML::Exception &cause) {
        throw New<Exception>(
          std::string("Failed to load file ") + fileName,
          SOURCE_LOCATION,
          New<Exception>(cause.what(), SourceLocation())
        );
      }
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
        throw New<Exception>(
          "unknown node type: " + yamlNode.Type(),
          SourceLocation(fileName, yamlNode.Mark().line)
        );
      }
    }

    Ptr<MapNode> parseMap(const YAML::Node &yamlNode) {
      auto node(New<MapNode>(SourceLocation(fileName, yamlNode.Mark().line)));
      for (auto iterator: yamlNode) {
        node->get(iterator.first.as<std::string>()) = parseNode(iterator.second);
      }
      return node;
    }

    Ptr<MapNode> parseSequence(const YAML::Node &yamlNode) {
      auto node(New<MapNode>(SourceLocation(fileName, yamlNode.Mark().line)));
      size_t index(0);
      for (auto subNode: yamlNode) {
        node->get(index++) = parseNode(subNode);
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
        tag.find(":Real") != std::string::npos
      ) {
        // default reals
        return parseAtom<Real<>>(yamlNode);
      } else if (tag.find(":Real64") != std::string::npos) {
        // explicit size reals
        return parseAtom<Real<64>>(yamlNode);
/*
      } else if (tag.find(":Real128") != std::string::npos) {
        // TODO: 128 bit real parsing
        // explicity size reals
        return parseAtom<Real<128>>(yamlNode);
*/
      } else if (tag.find(":Complex") != std::string::npos) {
        // default  complex
        return parseAtom<Complex<>>(yamlNode);
      } else if (tag.find(":Complex64") != std::string::npos) {
        // explicit size complex
        return parseAtom<Complex<64>>(yamlNode);
/*
      } else if (tag.find(":complex128") != std::string::npos) {
        // TODO: 128 bit real parsing
        // explicit size complex
        return parseAtom<Complex<128>>(yamlNode);
*/
      } else if (tag == "!" ){
        return parseAtom<std::string>(yamlNode);
      } else if (tag == "?" ){
        return parseImplicitTypeScalar(yamlNode);
      } else {
        auto value(yamlNode.Scalar());
        throw New<Exception>(
          "unsupported type tag " + tag + " for value " + value,
          SourceLocation(fileName, yamlNode.Mark().line)
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
          return New<SymbolNode>(
            yamlNode.Scalar(), SourceLocation(fileName, yamlNode.Mark().line)
          );
        }
      }
      // otherwise it's text
      return parseAtom<std::string>(yamlNode);
    }

    template <typename AtomicType>
    Ptr<AtomicNode<AtomicType>> parseAtom(const YAML::Node &yamlNode) {
      return New<AtomicNode<AtomicType>>(
        yamlNode.as<AtomicType>(),
        SourceLocation(fileName, yamlNode.Mark().line)
      );
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

