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

#include <Writer.hpp>
#include <Reader.hpp>
#include <util/Scanner.hpp>
#include <algorithms/Algorithm.hpp>
#include <tcc/Tcc.hpp>

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
      auto propertiesFileName(name + ".properties.yaml");
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
    static constexpr auto VERSION = "v1.0";

    template <typename F, typename TE>
    static Ptr<MapNode> writeTensor(
      const Ptr<Node> &node,
      const std::string &nodePath,
      const bool useBinary
    ) {
      auto pointerNode(node->toPtr<AtomicNode<Ptr<Object>>>());
      auto tensorExpression(
        dynamicPtrCast<TensorExpression<F,TE>>(pointerNode->value)
      );
      if (!tensorExpression) return nullptr;
      auto tensor(tensorExpression->evaluate());
      auto writtenTensor(New<MapNode>(SOURCE_LOCATION));

      // build dimensions from tensor length
      auto dimensions(New<MapNode>(SOURCE_LOCATION));
      for (size_t d(0); d < tensor->lens.size(); ++d) {
        auto dimension(New<MapNode>(SOURCE_LOCATION));
        dimension->setValue<>("length", tensor->lens[d]);
        if (tensor->dimensions[d]) {
          dimension->setValue<>("dimension", tensor->dimensions[d]->name);
        }
        dimensions->get(d) = dimension;
      }
      writtenTensor->get("dimensions") = dimensions;
      writtenTensor->setSymbol("scalarType", TypeTraits<F>::getName());
      writtenTensor->setValue<std::string>("version", VERSION);
      auto dataNode(New<MapNode>(SOURCE_LOCATION));
      dataNode->setValue<std::string>("type", "externalFile");
      dataNode->setValue<bool>("binary", useBinary);
      writtenTensor->get("data") = dataNode;

      // write tensor data as side effect
      if (useBinary) {
        writeTensorDataBinary(tensor, nodePath);
      } else {
        writeTensorDataText(tensor, nodePath);
      }

      return writtenTensor;
    }

    template <typename F, typename TE>
    static Ptr<AtomicNode<Ptr<Tensor<F,TE>>>> readTensor(
      const Ptr<MapNode> &node,
      const std::string &nodePath
    ) {
      auto version(node->getValue<std::string>("version"));
      ASSERT_LOCATION(
        version == VERSION,
        "Found and expected serialization versions mismatch. Found " +
          version + ", expecting " + VERSION,
        node->get("version")->sourceLocation
      );
      auto dataNode(node->getMap("data"));
      auto dataType(dataNode->getValue<std::string>("type"));
      ASSERT_LOCATION(
        dataType == "externalFile",
        "Only tensor data in external files are supported as 'data.type'",
        dataNode->get("type")->sourceLocation
      );
      auto dataBinary(dataNode->getValue<bool>("binary", false));
      auto fileName(nodePath + (dataBinary ? ".bin" : ".dat"));
      auto sourceLocation(node->sourceLocation);
      // get dimensions from meta data
      auto dimensions(node->getMap("dimensions"));
      std::vector<size_t> lens;
      for (auto key: dimensions->getKeys()) {
        lens.push_back(dimensions->getMap(key)->getValue<size_t>("length"));
      }
      if (dataBinary) {
        return readTensorDataBinary<F,TE>(fileName, lens, sourceLocation);
      } else {
        return readTensorDataText<F,TE>(fileName, lens, sourceLocation);
      }
    }

    template <typename F, typename TE>
    static void writeTensorDataText(
      const Ptr<Tensor<F,TE>> &tensor, const std::string &nodePath
    ) {
      constexpr size_t MAX_BUFFER_SIZE(128*1024*1024);
      std::string dataFileName(nodePath + ".dat");
      std::ofstream stream(dataFileName);
      if (stream.fail()) {
        std::stringstream explanation;
        explanation << "Failed to open file \"" << dataFileName << "\"";
        throw New<Exception>(explanation.str(), SOURCE_LOCATION);
      }
      // TODO: use stream Printer rather than <<
      stream << setprecision(17);

      // write the values only on root, all others still pariticipate calling MPI
      const size_t bufferSize(
        std::min(tensor->getElementsCount(), MAX_BUFFER_SIZE)
      );
      size_t localBufferSize(Cc4s::world->getRank() == 0 ? bufferSize : 0);
      std::vector<size_t> indices(localBufferSize);
      std::vector<F> values(localBufferSize);

      OUT() << "Writing to text file " << dataFileName << std::endl;
      if (Cc4s::options->dryRun) {
        stream << "# dry-run: no data written" << std::endl;
        return;
      }

      size_t index(0);
      LOG() << "indexCount=" << tensor->getElementsCount() << std::endl;
      while (index < tensor->getElementsCount()) {
        size_t elementsCount(
          std::min(bufferSize, tensor->getElementsCount()-index)
        );
        size_t localElementsCount(
          Cc4s::world->getRank() == 0 ? elementsCount : 0
        );
        for (size_t i(0); i < localElementsCount; ++i) {
          indices[i] = index+i;
        }
        LOG() << "reading " << elementsCount << " values from tensor..." << std::endl;
        tensor->read(localElementsCount, indices.data(), values.data());
        for (size_t i(0); i < localElementsCount; ++i) {
          stream << values[i] << "\n";
        }
        // wait until all processes finished
        Cc4s::world->barrier();
        index += elementsCount;
      }
      LOG() << "Written " << tensor->getElementsCount() <<
        " elements to text file " << dataFileName << std::endl;
    }

    template <typename F, typename TE>
    static void writeTensorDataBinary(
      const Ptr<Tensor<F,TE>> &tensor, const std::string &nodePath
    ) {
      std::string dataFileName(nodePath + ".bin");
      // open the file
      MPI_File file;
      int mpiError(
        MPI_File_open(
          Cc4s::world->getComm(), dataFileName.c_str(),
          MPI_MODE_CREATE | MPI_MODE_WRONLY,
          MPI_INFO_NULL, &file
        )
      );
      ASSERT_LOCATION(
        !mpiError, std::string("Failed to open file '") + dataFileName + "'",
        SOURCE_LOCATION
      )

      OUT() << "Writing to binary file " << dataFileName << std::endl;
      if (Cc4s::options->dryRun) return;

      // write tensor data with values from file
      tensor->readToFile(file);
      LOG() << "Written " << sizeof(F)*tensor->getElementsCount() <<
        " bytes to binary file " << dataFileName << std::endl;

      // done
      MPI_File_close(&file);
    }

    template <typename F, typename TE>
    static Ptr<AtomicNode<Ptr<Tensor<F,TE>>>> readTensorDataText(
      const std::string &fileName,
      const std::vector<size_t> &lens,
      const SourceLocation &sourceLocation
    ) {
      constexpr size_t MAX_BUFFER_SIZE(128*1024*1024);
      std::ifstream stream(fileName.c_str());
      if (stream.fail()) {
        std::stringstream explanation;
        explanation << "Failed to open file \"" << fileName << "\"";
        throw New<Exception>(explanation.str(), sourceLocation);
      }
      Scanner scanner(&stream);

      // create tensor
      auto A( Tcc<TE>::template tensor<F>(lens, fileName) );
      auto result(
        New<AtomicNode<Ptr<Tensor<F,TE>>>>(A, SourceLocation(fileName,1))
      );
      OUT() << "Reading from text file " << fileName << std::endl;
      if (Cc4s::options->dryRun) return result;

      // read the values only on root, all others still pariticipate calling MPI
      const size_t bufferSize(std::min(A->getElementsCount(), MAX_BUFFER_SIZE));
      size_t localBufferSize(Cc4s::world->getRank() == 0 ? bufferSize : 0);
      std::vector<size_t> indices(localBufferSize);
      std::vector<F> values(localBufferSize);

      size_t index(0);
      LOG() << "indexCount=" << A->getElementsCount() << std::endl;
      NumberScanner<F> numberScanner(&scanner);
      while (index < A->getElementsCount()) {
        size_t elementsCount(std::min(bufferSize, A->getElementsCount()-index));
        size_t localElementsCount(Cc4s::world->getRank() == 0 ? elementsCount : 0);
        for (size_t i(0); i < localElementsCount; ++i) {
          indices[i] = index+i;
          values[i] = numberScanner.nextNumber();
        }
        // wait until all processes finished reading this buffer into the tensor
        Cc4s::world->barrier();
        LOG() << "writing " << elementsCount << " values to tensor..." << std::endl;
        A->write(localElementsCount, indices.data(), values.data());
        index += elementsCount;
      }
      LOG() << "Read " << A->getElementsCount() <<
        " elements from text file " << fileName << std::endl;

      return result;
    }

    template <typename F, typename TE>
    static Ptr<AtomicNode<Ptr<Tensor<F,TE>>>> readTensorDataBinary(
      const std::string &fileName,
      const std::vector<size_t> &lens,
      const SourceLocation &sourceLocation
    ) {
      // open the file
      MPI_File file;
      int mpiError(
        MPI_File_open(
          Cc4s::world->getComm(), fileName.c_str(), MPI_MODE_RDONLY,
          MPI_INFO_NULL, &file
        )
      );
      ASSERT_LOCATION(
        !mpiError, std::string("Failed to open file '") + fileName + "'",
        sourceLocation
      )

      // create tensor
      auto A( Tcc<TE>::template tensor<F>(lens, fileName) );
      auto result(
        New<AtomicNode<Ptr<Tensor<F,TE>>>>(A, SourceLocation(fileName,1))
      );

      OUT() << "Reading from binary file " << fileName << std::endl;
      if (Cc4s::options->dryRun) return result;
      // write tensor data with values from file
      A->writeFromFile(file);
      LOG() << "Read " << sizeof(F)*A->getElementsCount() <<
        " bytes from binary file " << fileName << std::endl;

      // done
      MPI_File_close(&file);
      return result;
    }

    static int WRITE_REGISTERED, READ_REGISTERED;
  };
}

#endif

