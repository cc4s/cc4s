/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
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

ALGORITHM_REGISTRAR_DEFINITION(TensorReader);

/**
 * \brief Reader for tensor data.
 */
Ptr<MapNode> TensorReader::run(const Ptr<MapNode> &arguments) {
  auto fileName(arguments->getValue<std::string>("fileName"));
  // get tensor meta data from given file
  auto tensor(Parser(fileName).parse()->map());
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
  auto dimensions(tensor->getMap("dimensions"));
  std::vector<size_t> lens;
  for (auto key: dimensions->getKeys()) {
    lens.push_back(dimensions->getMap(key)->getValue<size_t>("length"));
  }

  Ptr<Node> tensorData;
  // multiplex calls to template methods depending on tensor engine and type
  if (binary) {
    if (Cc4s::options->dryRun) {
      if (scalarType == "real64") {
        tensorData = readBinary<Real<64>,DryTensorEngine>(
          fileName, lens, sourceLocation
        );
      } else if (scalarType == "complex64") {
        tensorData = readBinary<Complex<64>,DryTensorEngine>(
          fileName, lens, sourceLocation
        );
      } else {
        ASSERT_LOCATION(
          false, "scalar type '" + scalarType + "' not supported", sourceLocation
        );
      }
    } else {
      if (scalarType == "real64") {
        tensorData = readBinary<Real<64>,DefaultTensorEngine>(
          fileName, lens, sourceLocation
        );
      } else if (scalarType == "complex64") {
        tensorData = readBinary<Complex<64>,DefaultTensorEngine>(
          fileName, lens, sourceLocation
        );
      } else {
        ASSERT_LOCATION(
          false, "scalar type '" + scalarType + "' not supported", sourceLocation
        );
      }
    }
  } else {
    if (Cc4s::options->dryRun) {
      if (scalarType == "real64") {
        tensorData = readText<Real<64>,DryTensorEngine>(
          fileName, lens, sourceLocation
        );
      } else if (scalarType == "complex64") {
        tensorData = readText<Complex<64>,DryTensorEngine>(
          fileName, lens, sourceLocation
        );
      } else {
        ASSERT_LOCATION(
          false, "scalar type '" + scalarType + "' not supported", sourceLocation
        );
      }
    } else {
      if (scalarType == "real64") {
        tensorData = readText<Real<64>,DefaultTensorEngine>(
          fileName, lens, sourceLocation
        );
      } else if (scalarType == "complex64") {
        tensorData = readText<Complex<64>,DefaultTensorEngine>(
          fileName, lens, sourceLocation
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

  return New<AtomicNode<Ptr<Tensor<F,TE>>>>(A, SourceLocation(fileName,1));
}

template <typename F, typename TE>
Ptr<AtomicNode<Ptr<Tensor<F,TE>>>> TensorReader::readBinary(
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

  // write tensor data with values from file
  A->writeFromFile(file);

  // done
  MPI_File_close(&file);
  return New<AtomicNode<Ptr<Tensor<F,TE>>>>(A, SourceLocation(fileName,1));
}

