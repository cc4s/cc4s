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

#ifndef TENSOR_IO_DEFINED
#define TENSOR_IO_DEFINED

#include <algorithms/Algorithm.hpp>
#include <tcc/Tcc.hpp>
#include <SharedPointer.hpp>
#include <Writer.hpp>
#include <Reader.hpp>
#include <Scanner.hpp>

namespace cc4s {
  class TensorIo {
  public:
    // TODO: use context to hold dimension info
    static Ptr<TensorDimension> getDimension(const std::string &name) {
      // check if name is entetered in map
      auto iterator(TensorDimension::dimensions.find(name));
      if (iterator != TensorDimension::dimensions.end()) return iterator->second;
      // otherwise: create new tensor dimension entry
      auto tensorDimension(New<TensorDimension>());
      tensorDimension->name = name;
      TensorDimension::dimensions[name] = tensorDimension;

      // check if file exists specifying dimension properties
      auto propertiesFileName(name + ".yaml");
      if (std::ifstream(propertiesFileName).good()) {
        auto propertiesMap(
          Parser(propertiesFileName).parse()->toPtr<MapNode>()
        );
        for (auto key: propertiesMap->getKeys()) {
          auto propertyMap(propertiesMap->getMap(key));
          auto property(New<TensorDimensionProperty>());
          property->name = key;
          for (auto indexKey: propertyMap->getKeys()) {
            auto index(std::stol(indexKey));
            auto propertyValue(propertyMap->getValue<Natural<>>(indexKey));
            // enter index -> property map
            property->propertyOfIndex[index] = propertyValue;
            // build reverse lookup map of sets as well
            property->indicesOfProperty[propertyValue].insert(index);
          }
          LOG() << "entering property "
            << property->name << " of dimension " << name << std::endl;
          tensorDimension->properties[property->name] = property;
        }
      }
      return tensorDimension;
    }

    /**
     * \brief Static handler routine for writing nodes of the type
     * AtomicNode<Ptr<Tensor<F,TE>>>. If the given node is of such a type
     * the contained tensor will be written to disk and a MapNode is
     * constructed and returned containing all necessary information to
     * load the written tensor again.
     * nullptr is returned if the given node is not an AtomicNode containing
     * a tensor pointer.
     **/
    static Ptr<Node> write(
      const Ptr<Node> &node, const std::string &nodePath, const bool useBinary
    ) {
      // multiplex different tensor types
      Ptr<Node> writtenNode;
      if (!Cc4s::options->dryRun) {
        using TE = DefaultTensorEngine;
        writtenNode = writeTensor<Real<64>,TE>(node, nodePath, useBinary);
        if (writtenNode) return writtenNode;
        writtenNode = writeTensor<Complex<64>,TE>(node, nodePath, useBinary);
        if (writtenNode) return writtenNode;
/* // TODO: fully support 128 bit floats or remove
        writtenNode = writeTensor<Real<128>,TE>(node, nodePath);
        if (writtenNode) return writtenNode;
        writtenNode = writeTensor<Complex<128>,TE>(node, nodePath);
        if (writtenNode) return writtenNode;
*/
      } else {
        using TE = DefaultDryTensorEngine;
        writtenNode = writeTensor<Real<64>,TE>(node, nodePath, useBinary);
        if (writtenNode) return writtenNode;
        writtenNode = writeTensor<Complex<64>,TE>(node, nodePath, useBinary);
        if (writtenNode) return writtenNode;
/*
        writtenNode = writeTensor<Real<128>,TE>(node, nodePath);
        if (writtenNode) return writtenNode;
        writtenNode = writeTensor<Complex<128>,TE>(node, nodePath);
        if (writtenNode) return writtenNode;
*/
      }
      // otherwise, not my type ... let another write routine handle this data
      return nullptr;
    }

    static Ptr<Node> read(
      const Ptr<MapNode> &node, const std::string &nodePath
    ) {
      auto scalarType(node->getValue<std::string>("scalarType"));
      // multiplex different tensor types
      Ptr<Node> workingNode;
      if (!Cc4s::options->dryRun) {
        using TE = DefaultTensorEngine;
        if (scalarType == TypeTraits<Real<64>>::getName()) {
          return readTensor<Real<64>,TE>(node, nodePath);
        } else if (scalarType == TypeTraits<Complex<64>>::getName()) {
          return readTensor<Complex<64>,TE>(node, nodePath);
        }
      } else {
        using TE = DefaultDryTensorEngine;
        if (scalarType == TypeTraits<Real<64>>::getName()) {
          return readTensor<Real<64>,TE>(node, nodePath);
        } else if (scalarType == TypeTraits<Complex<64>>::getName()) {
          return readTensor<Complex<64>,TE>(node, nodePath);
        }
      }
      std::stringstream explanation;
      explanation << "Unsupported scalarType \"" << scalarType << "\"";
      throw New<Exception>(explanation.str(), SOURCE_LOCATION);
    }

  protected:
    /**
     * \brief Serialization version. Only objects written by
     * a matching serialization version can be read. Increase this version
     * if the data written is no longer backward compatible.
     **/
    static constexpr auto VERSION = "1.0";

    template <typename F, typename TE>
    static Ptr<MapNode> writeTensor(
      const Ptr<Node> &node,
      const std::string &nodePath,
      const bool useBinary
    );

    template <typename F, typename TE>
    static Ptr<AtomicNode<Ptr<Tensor<F,TE>>>> readTensor(
      const Ptr<MapNode> &node,
      const std::string &nodePath
    );

    template <typename F, typename TE>
    static void writeTensorElementsText(
      const Ptr<Tensor<F,TE>> &tensor, const std::string &nodePath
    );

    template <typename F, typename TE>
    static void writeTensorElementsBinary(
      const Ptr<Tensor<F,TE>> &tensor, const std::string &nodePath
    );

    template <typename F, typename TE>
    static Ptr<AtomicNode<Ptr<Tensor<F,TE>>>> readTensorElementsText(
      const std::string &fileName,
      const std::vector<size_t> &lens,
      const SourceLocation &sourceLocation
    );

    template <typename F, typename TE>
    static Ptr<AtomicNode<Ptr<Tensor<F,TE>>>> readTensorElementsBinary(
      const std::string &fileName,
      const std::vector<size_t> &lens,
      const SourceLocation &sourceLocation
    );

    static int WRITE_REGISTERED, READ_REGISTERED;
  };
}

#endif

