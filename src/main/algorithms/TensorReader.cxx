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

#include <algorithms/TensorReader.hpp>

#include <tcc/TensorExpression.hpp>
#include <Node.hpp>
#include <Cc4s.hpp>
#include <Parser.hpp>
#include <SharedPointer.hpp>
#include <Scanner.hpp>
#include <TensorIo.hpp>
#include <Real.hpp>
#include <Complex.hpp>
#include <MathFunctions.hpp>

#include <vector>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(TensorReader)

/**
 * \brief Reader for tensor data.
 */
Ptr<MapNode> TensorReader::run(const Ptr<MapNode> &arguments) {
  auto fileName(arguments->getValue<std::string>("fileName"));
  // get tensor meta data from given file
  auto tensor(Parser(fileName).parse()->toPtr<MapNode>());
  auto scalarType(tensor->getValue<std::string>("scalarType"));
  auto binary(tensor->getValue<bool>("binary", false));

  // if fileName contains '/' change directory
  char currentDirectory[PATH_MAX];
  getcwd(currentDirectory, sizeof(currentDirectory));
  // TODO: support other filesystems
  auto slashPosition(fileName.rfind('/'));
  if (slashPosition != std::string::npos) {
    auto fileDirectory(fileName.substr(0, slashPosition));
    chdir(fileDirectory.c_str());
  }

  // read tensor data from file named by 'data' entry
  if (tensor->get("components")) {
    auto components(tensor->getMap("components"));
    for (auto key: components->getKeys()) {
      readData(components->getMap(key), scalarType, binary);
    }
  } else {
    readData(tensor, scalarType, binary);
  }

  // change back to current directory
  chdir(currentDirectory);

  // create result
  auto result(New<MapNode>(tensor->sourceLocation));
  result->get("tensor") = tensor;
  return result;
}

void TensorReader::readData(
  const Ptr<MapNode> &tensor,
  const std::string &scalarType,
  const bool binary
) {
  auto fileName(tensor->getValue<std::string>("data"));
  auto sourceLocation(tensor->sourceLocation);
  // get dimensions from meta data
  auto dimensionsMap(tensor->getMap("dimensions"));
  std::vector<size_t> lens;
  std::vector<Ptr<TensorDimension>> dimensions;
  for (auto key: dimensionsMap->getKeys()) {
    lens.push_back(dimensionsMap->getMap(key)->getValue<Natural<>>("length"));
    // get dimension properties, if given
    dimensions.push_back(
      TensorIo::getDimension(
        dimensionsMap->getMap(key)->getValue<std::string>("type")
      )
    );
  }

  // check if tensor is given in block sparse format
  auto nonZeroConditions(New<TensorNonZeroConditions>());
  auto nonZeroConditionsNode(tensor->get("nonZeroCondition"));
  if (nonZeroConditionsNode) {
    OUT() << "Unpacking block-sparse tensor" << std::endl;
    // currently only "all" is supported as non-zero conditions
    auto allMap(nonZeroConditionsNode->toPtr<MapNode>()->getMap("all"));
    for (auto key: allMap->getKeys()) {
      auto nonZeroCondition(New<TensorNonZeroCondition>());
      nonZeroCondition->name = allMap->getMap(key)->getValue<std::string>("name");
      auto propertiesMap(allMap->getMap(key)->getMap("properties"));
      Natural<> D(0);
      for (auto propertyKey: propertiesMap->getKeys()) {
        auto propertyMap(propertiesMap->getMap(propertyKey));
        TensorDimensionPropertyReference dimensionProperty;
        dimensionProperty.dimension = propertyMap->getValue<Natural<>>("dimension");
        auto dimension(dimensions[dimensionProperty.dimension]);
        dimensionProperty.property = dimension->properties[
          propertyMap->getValue<std::string>("property")
        ];
        nonZeroCondition->dimensionPropertyReferences.push_back(dimensionProperty);
        ++D;
      }
      auto tuplesMap(allMap->getMap(key)->getMap("nonZeros"));
      for (Natural<> i(0); i < tuplesMap->size(); ++i) {
        std::vector<Natural<>> tuple(D);
        for (Natural<> d(0); d < D; ++d) {
          tuple[d] = tuplesMap->getMap(i)->getValue<Natural<>>(d);
        }
        nonZeroCondition->tuples.push_back(tuple);
      }
      nonZeroConditions->all.push_back(nonZeroCondition);
    }
  }

  Ptr<Node> tensorData;
  // multiplex calls to template methods depending on tensor engine and type
  if (binary) {
    if (Cc4s::options->dryRun) {
      using TE = DefaultDryTensorEngine;
      if (scalarType == "real64") {
        tensorData = readBinary<Real<64>,TE>(
          fileName, lens, dimensions, nonZeroConditions, sourceLocation
        );
      } else if (scalarType == "complex64") {
        tensorData = readBinary<Complex<64>,TE>(
          fileName, lens, dimensions, nonZeroConditions, sourceLocation
        );
      } else {
        ASSERT_LOCATION(
          false, "scalar type '" + scalarType + "' not supported", sourceLocation
        );
      }
    } else {
      using TE = DefaultTensorEngine;
      if (scalarType == "real64") {
        tensorData = readBinary<Real<64>,TE>(
          fileName, lens, dimensions, nonZeroConditions, sourceLocation
        );
      } else if (scalarType == "complex64") {
        tensorData = readBinary<Complex<64>,TE>(
          fileName, lens, dimensions, nonZeroConditions, sourceLocation
        );
      } else {
        ASSERT_LOCATION(
          false, "scalar type '" + scalarType + "' not supported", sourceLocation
        );
      }
    }
  } else {
    if (Cc4s::options->dryRun) {
      using TE = DefaultDryTensorEngine;
      if (scalarType == "real64") {
        tensorData = readText<Real<64>,TE>(
          fileName, lens, dimensions, nonZeroConditions, sourceLocation
        );
      } else if (scalarType == "complex64") {
        tensorData = readText<Complex<64>,TE>(
          fileName, lens, dimensions, nonZeroConditions, sourceLocation
        );
      } else {
        ASSERT_LOCATION(
          false, "scalar type '" + scalarType + "' not supported", sourceLocation
        );
      }
    } else {
      using TE = DefaultTensorEngine;
      if (scalarType == "real64") {
        tensorData = readText<Real<64>,TE>(
          fileName, lens, dimensions, nonZeroConditions, sourceLocation
        );
      } else if (scalarType == "complex64") {
        tensorData = readText<Complex<64>,TE>(
          fileName, lens, dimensions, nonZeroConditions, sourceLocation
        );
      } else {
        ASSERT_LOCATION(
          false, "scalar type '" + scalarType + "' not supported", sourceLocation
        );
      }
    }
  }

  // replace data entry with actual tensor
  tensor->get("data") = tensorData;
}

template <typename F, typename TE>
Ptr<Node> TensorReader::readText(
  const std::string &fileName,
  const std::vector<size_t> &lens,
  const std::vector<Ptr<TensorDimension>> &dimensions,
  const Ptr<TensorNonZeroConditions> &nonZeroConditions,
  const SourceLocation &sourceLocation
) {
  constexpr size_t MAX_BUFFER_SIZE(32*1024*1024);
  std::ifstream stream(fileName.c_str());
  if (stream.fail()) {
    std::stringstream explanation;
    explanation << "Failed to open file \"" << fileName << "\"";
    throw New<Exception>(explanation.str(), sourceLocation);
  }
  Scanner scanner(&stream);

  // create dense result tensor
  auto A( Tcc<TE>::template tensor<F>(lens, fileName) );
  // TODO: API for accessing tensor dimension info
  A->dimensions = dimensions;
  A->nonZeroConditions = nonZeroConditions;
  auto result(
    New<PointerNode<TensorExpression<F,TE>>>(A, SourceLocation(fileName,1))
  );
  OUT() << "Reading from text file " << fileName << std::endl;
  if (Cc4s::options->dryRun) return result;

  LOG() << "#non-zero conditions: " << nonZeroConditions->all.size() << std::endl;
  TensorNonZeroBlockIterator blockIterator(nonZeroConditions);
  while (!blockIterator.atEnd()) {
    auto elementIterator(blockIterator.getElementIterator(lens));
    auto blockElementsCount(elementIterator.getCount());
    LOG() << "reading block:" << blockIterator.print() <<
      " with " << blockElementsCount << " elements" << std::endl;

    // read the values only on root, all others still pariticipate calling MPI
    const size_t bufferSize(std::min(blockElementsCount, MAX_BUFFER_SIZE));
    size_t localBufferSize(Cc4s::world->getRank() == 0 ? bufferSize : 0);
    std::vector<size_t> indices(localBufferSize);
    std::vector<F> values(localBufferSize);

    size_t readElementsCount(0);
    NumberScanner<F> numberScanner(&scanner);
    while (readElementsCount < blockElementsCount) {
      size_t elementsCount(
        std::min(bufferSize, blockElementsCount - readElementsCount)
      );
      size_t localElementsCount(Cc4s::world->getRank() == 0 ? elementsCount : 0);
      for (size_t i(0); i < localElementsCount; ++i) {
        indices[i] = elementIterator.getGlobalIndex();
        values[i] = numberScanner.nextNumber();
        ++elementIterator;
      }
      // wait until all processes finished reading this buffer into the tensor
      Cc4s::world->barrier();
      LOG() << "writing " << elementsCount << " values to tensor..." << std::endl;
      A->write(localElementsCount, indices.data(), values.data());
      readElementsCount += elementsCount;
    }

    ++blockIterator;
  }

  return result;
}

template <typename F, typename TE>
Ptr<Node> TensorReader::readBinary(
  const std::string &fileName,
  const std::vector<size_t> &lens,
  const std::vector<Ptr<TensorDimension>> &dimensions,
  const Ptr<TensorNonZeroConditions> &nonZeroConditions,
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
  A->nonZeroConditions = nonZeroConditions;
  auto result(
    New<PointerNode<TensorExpression<F,TE>>>(A, SourceLocation(fileName,1))
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

