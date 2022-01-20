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

#include <TensorIo.hpp>
#include <Reader.hpp>

using namespace cc4s;

bool TensorIo::WRITE_REGISTERED =
  Writer::registerWriteFunction("Tensor", TensorIo::write);

bool TensorIo::READ_REGISTERED =
  Reader::registerReadFunction("Tensor", TensorIo::read);

const Natural<32> TensorIo::VERSION = 100;


Ptr<TensorDimension> TensorIo::getDimension(const std::string &name) {
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

Ptr<Node> TensorIo::write(
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
  } else {
    using TE = DefaultDryTensorEngine;
    writtenNode = writeTensor<Real<64>,TE>(node, nodePath, useBinary);
    if (writtenNode) return writtenNode;
    writtenNode = writeTensor<Complex<64>,TE>(node, nodePath, useBinary);
    if (writtenNode) return writtenNode;
  }
  // otherwise, not my type ... let another write routine handle this data
  return nullptr;
}

Ptr<Node> TensorIo::read(
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


template <typename F, typename TE>
Ptr<MapNode> TensorIo::writeTensor(
  const Ptr<Node> &node,
  const std::string &nodePath,
  const bool useBinary
) {
  auto pointerNode(node->toPtr<AtomicNode<Ptr<Object>>>());
  if (!pointerNode) return nullptr;
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
    dimension->setValue("length", tensor->lens[d]);
    if (tensor->dimensions.size() > 0 && tensor->dimensions[d]) {
      dimension->setValue("type", tensor->dimensions[d]->name);
    }
    dimensions->get(d) = dimension;
  }
  writtenTensor->get("dimensions") = dimensions;
  writtenTensor->setSymbol("scalarType", TypeTraits<F>::getName());
  writtenTensor->setValue("version", VERSION);
  auto elementsNode(New<MapNode>(SOURCE_LOCATION));
  writtenTensor->get("elements") = elementsNode;
  writtenTensor->setValue("unit", tensor->getUnit());
  writtenTensor->get("metaData") = tensor->getMetaData();

  // write tensor elements as side effect
  if (useBinary) {
    elementsNode->setValue("type", std::string("IeeeBinaryFile"));
    writeTensorElementsBinary(tensor, nodePath);
  } else {
    elementsNode->setValue("type", std::string("TextFile"));
    writeTensorElementsText(tensor, nodePath);
  }

  return writtenTensor;
}

template <typename F, typename TE>
Ptr<PointerNode<Tensor<F,TE>>> TensorIo::readTensor(
  const Ptr<MapNode> &node,
  const std::string &nodePath
) {
  auto version(node->getValue<Natural<32>>("version"));
  ASSERT_LOCATION(
    version == VERSION,
    "Found and expected serialization versions mismatch. Found " +
      std::to_string(version) + ", expecting " + std::to_string(VERSION),
    node->get("version")->sourceLocation
  );
  auto elementsNode(node->getMap("elements"));
  auto elementsType(elementsNode->getValue<std::string>("type"));
  auto elementsBinary(elementsType == "IeeBinaryFile" ? true : false);
  auto elementsPath(nodePath + ".elements");
  auto sourceLocation(node->sourceLocation);
  auto dimensions(node->getMap("dimensions"));

  std::vector<size_t> lens;
  for (auto key: dimensions->getKeys()) {
    lens.push_back(dimensions->getMap(key)->getValue<size_t>("length"));
  }
  Ptr<Tensor<F,TE>> tensor;
  if (elementsBinary) {
    tensor = readTensorElementsBinary<F,TE>(elementsPath, lens, sourceLocation);
  } else {
    tensor = readTensorElementsText<F,TE>(elementsPath, lens, sourceLocation);
  }
  tensor->unit = node->getValue<Real<>>("unit");
  tensor->metaData = node->getMap("metaData");

  return New<PointerNode<Tensor<F,TE>>>(tensor, SourceLocation(nodePath,1));
}

template <typename F, typename TE>
void TensorIo::writeTensorElementsText(
  const Ptr<Tensor<F,TE>> &tensor, const std::string &nodePath
) {
  // TODO: query max buffer size for tensor write from tensor engine
  constexpr size_t MAX_BUFFER_SIZE(32*1024*1024);
  std::string elementsPath(nodePath + ".elements");
  std::ofstream stream(elementsPath);
  if (stream.fail()) {
    std::stringstream explanation;
    explanation << "Failed to open file \"" << elementsPath << "\"";
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

  OUT() << "Writing to text file " << elementsPath << std::endl;
  if (Cc4s::options->dryRun) {
    stream << "# dry-run: no elements written" << std::endl;
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
    " elements to text file " << elementsPath << std::endl;
}

template <typename F, typename TE>
void TensorIo::writeTensorElementsBinary(
  const Ptr<Tensor<F,TE>> &tensor, const std::string &nodePath
) {
  std::string elementsPath(nodePath + ".elements");
  // open the file
  MPI_File file;
  int mpiError(
    MPI_File_open(
      Cc4s::world->getComm(), elementsPath.c_str(),
      MPI_MODE_CREATE | MPI_MODE_WRONLY,
      MPI_INFO_NULL, &file
    )
  );
  ASSERT_LOCATION(
    !mpiError, std::string("Failed to open file '") + elementsPath + "'",
    SOURCE_LOCATION
  )

  OUT() << "Writing to binary file " << elementsPath << std::endl;
  if (Cc4s::options->dryRun) return;

  // write tensor elements with values from file
  tensor->readToFile(file);
  LOG() << "Written " << sizeof(F)*tensor->getElementsCount() <<
    " bytes to binary file " << elementsPath << std::endl;

  // done
  MPI_File_close(&file);
}

template <typename F, typename TE>
Ptr<Tensor<F,TE>> TensorIo::readTensorElementsText(
  const std::string &fileName,
  const std::vector<size_t> &lens,
  const SourceLocation &sourceLocation
) {
  // TODO: query max buffer size for tensor write from tensor engine
  constexpr size_t MAX_BUFFER_SIZE(32*1024*1024);
  std::ifstream stream(fileName.c_str());
  if (stream.fail()) {
    std::stringstream explanation;
    explanation << "Failed to open file \"" << fileName << "\"";
    throw New<Exception>(explanation.str(), sourceLocation);
  }
  Scanner scanner(&stream);

  // create tensor
  auto tensor( Tcc<TE>::template tensor<F>(lens, fileName) );
  OUT() << "Reading from text file " << fileName << std::endl;
  if (Cc4s::options->dryRun) return tensor;

  // read the values only on root, all others still pariticipate calling MPI
  const size_t bufferSize(std::min(tensor->getElementsCount(), MAX_BUFFER_SIZE));
  size_t localBufferSize(Cc4s::world->getRank() == 0 ? bufferSize : 0);
  std::vector<size_t> indices(localBufferSize);
  std::vector<F> values(localBufferSize);

  size_t index(0);
  LOG() << "indexCount=" << tensor->getElementsCount() << std::endl;
  NumberScanner<F> numberScanner(&scanner);
  while (index < tensor->getElementsCount()) {
    size_t elementsCount(std::min(bufferSize, tensor->getElementsCount()-index));
    size_t localElementsCount(Cc4s::world->getRank() == 0 ? elementsCount : 0);
    for (size_t i(0); i < localElementsCount; ++i) {
      indices[i] = index+i;
      values[i] = numberScanner.nextNumber();
    }
    // wait until all processes finished reading this buffer into the tensor
    Cc4s::world->barrier();
    LOG() << "writing " << elementsCount << " values to tensor..." << std::endl;
    tensor->write(localElementsCount, indices.data(), values.data());
    index += elementsCount;
  }
  LOG() << "Read " << tensor->getElementsCount() <<
    " elements from text file " << fileName << std::endl;

  return tensor;
}

template <typename F, typename TE>
Ptr<Tensor<F,TE>> TensorIo::readTensorElementsBinary(
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
  auto tensor( Tcc<TE>::template tensor<F>(lens, fileName) );

  OUT() << "Reading from binary file " << fileName << std::endl;
  if (Cc4s::options->dryRun) return tensor;
  // write tensor elements with values from file
  tensor->writeFromFile(file);
  LOG() << "Read " << sizeof(F)*tensor->getElementsCount() <<
    " bytes from binary file " << fileName << std::endl;

  // done
  MPI_File_close(&file);
  return tensor;
}

