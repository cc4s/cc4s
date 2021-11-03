#include <algorithms/TensorReader.hpp>

#include <Cc4s.hpp>
#include <Parser.hpp>
#include <util/Scanner.hpp>
#include <tcc/Tcc.hpp>
#include <math/Real.hpp>
#include <math/Complex.hpp>
#include <math/MathFunctions.hpp>
#include <vector>

using namespace cc4s;

// TODO: table of tensor dimension information, belongs in tcc
std::map<std::string, Ptr<TensorDimension>> tensorDimensions;

ALGORITHM_REGISTRAR_DEFINITION(TensorReader)

/**
 * \brief Reader for tensor data.
 */
Ptr<MapNode> TensorReader::run(const Ptr<MapNode> &arguments) {
  auto fileName(arguments->getValue<std::string>("fileName"));
  // get tensor meta data from given file
  auto tensor(Parser(fileName).parse()->toMap());
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

// TODO: belongs in tcc
Ptr<TensorDimension> getDimension(const std::string &name) {
  // check if name is entetered in map
  auto iterator(tensorDimensions.find(name));
  if (iterator != tensorDimensions.end()) return iterator->second;
  // otherwise: create new tensor dimension entry
  auto tensorDimension(New<TensorDimension>());
  tensorDimension->name = name;
  // check if file exists specifying dimension properties
  try {
    auto propertiesMap(Parser(name + ".properties.yaml").parse()->toMap());
    for (auto key: propertiesMap->getKeys()) {
      auto propertyMap(propertiesMap->getMap(key));
      auto property(New<TensorDimensionProperty>());
      property->name = key;
      for (auto indexKey: propertyMap->getKeys()) {
        auto index(std::stol(indexKey));
        auto propertyValue(propertyMap->getValue<Natural>(indexKey));
        // enter index -> property map
        property->propertyOfIndex[index] = propertyValue;
        // build reverse lookup map of sets as well
        property->indicesOfProperty[propertyValue].insert(index);
      }
      LOG() << "entering property " << property->name << " of dimension " << name << std::endl;
      tensorDimension->properties[property->name] = property;
    }
  } catch (Ptr<Exception> &cause) {
    throw cause;
    // no properties file could be loaded: dimension without properties
  }
  return tensorDimensions[name] = tensorDimension;
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
    lens.push_back(dimensionsMap->getMap(key)->getValue<size_t>("length"));
    // get dimension properties, if given
    dimensions.push_back(
      getDimension(dimensionsMap->getMap(key)->getValue<std::string>("type"))
    );
  }

  // check if tensor is given in block sparse format
  auto nonZeroConditions(New<TensorNonZeroConditions>());
  auto nonZeroConditionsNode(tensor->get("nonZeroCondition"));
  if (nonZeroConditionsNode) {
    OUT() << "Unpacking block-sparse tensor" << std::endl;
    // currently only "all" is supported as non-zero conditions
    auto allMap(nonZeroConditionsNode->toMap()->getMap("all"));
    for (auto key: allMap->getKeys()) {
      auto nonZeroCondition(New<TensorNonZeroCondition>());
      auto propertiesMap(allMap->getMap(key)->getMap("properties"));
      Natural D(0);
      for (auto propertyKey: propertiesMap->getKeys()) {
        auto propertyMap(propertiesMap->getMap(propertyKey));
        TensorDimensionPropertyReference dimensionProperty;
        dimensionProperty.dimension = propertyMap->getValue<Natural>("dimension");
        auto dimension(dimensions[dimensionProperty.dimension]);
        dimensionProperty.property = dimension->properties[
          propertyMap->getValue<std::string>("property")
        ];
        nonZeroCondition->dimensionPropertyReferences.push_back(dimensionProperty);
        ++D;
      }
      auto tuplesMap(allMap->getMap(key)->getMap("nonZeros"));
      for (Natural i(0); i < tuplesMap->size(); ++i) {
        std::vector<Natural> tuple(D);
        for (Natural d(0); d < D; ++d) {
          tuple[d] = tuplesMap->getMap(i)->getValue<Natural>(d);
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
Ptr<AtomicNode<Ptr<Tensor<F,TE>>>> TensorReader::readText(
  const std::string &fileName,
  const std::vector<size_t> &lens,
  const std::vector<Ptr<TensorDimension>> &dimensions,
  const Ptr<TensorNonZeroConditions> &nonZeroConditions,
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

  // create dense result tensor
  auto A( Tcc<TE>::template tensor<F>(lens, fileName) );
  auto result(
    New<AtomicNode<Ptr<Tensor<F,TE>>>>(A, SourceLocation(fileName,1))
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
Ptr<AtomicNode<Ptr<Tensor<F,TE>>>> TensorReader::readBinary(
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

