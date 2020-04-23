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
  // get dimensions from meta data
  auto dimensions(tensor->getMap("dimensions"));
  std::vector<size_t> lens;
  for (auto key: dimensions->getKeys()) {
    lens.push_back(dimensions->getValue<size_t>(key));
  }

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
  auto tensorData(
    readText(
      tensor->getValue<std::string>("data"),
      lens,
      tensor->getValue<std::string>("scalarType")
    )
  );

  // change back to current directory
  chdir(currentDirectory);

  // replace data entry with actual tensor
  tensor->get("data") = tensorData;

  // create result
  auto result(New<MapNode>());
  result->get("tensor") = tensor;
  return result;
}

Ptr<Node> TensorReader::readText(
  const std::string &fileName,
  const std::vector<size_t> &lens,
  const std::string &scalarType
) {
  // multiplex calls to template methods depending on tensor engine and type
  if (Cc4s::options->dryRun) {
    if (scalarType == "real") {
      return readText<Real<>,DryTensorEngine>(fileName, lens);
    } else if (scalarType == "real64") {
      return readText<Real<64>,DryTensorEngine>(fileName, lens);
    } else if (scalarType == "complex") {
      return readText<Complex<>,DryTensorEngine>(fileName, lens);
    } else if (scalarType == "complex64") {
      return readText<Complex<64>,DryTensorEngine>(fileName, lens);
    } else {
      Assert(false, "scalar type '" + scalarType + "' not supported");
    }
  } else {
    if (scalarType == "real") {
      return readText<Real<>,DefaultTensorEngine>(fileName, lens);
    } else if (scalarType == "real64") {
      return readText<Real<64>,DefaultTensorEngine>(fileName, lens);
    } else if (scalarType == "complex") {
      return readText<Complex<>,DefaultTensorEngine>(fileName, lens);
    } else if (scalarType == "complex64") {
      return readText<Complex<64>,DefaultTensorEngine>(fileName, lens);
    } else {
      Assert(false, "scalar type '" + scalarType + "' not supported");
    }
  }
}

template <typename F, typename TE>
Ptr<AtomicNode<Ptr<Tensor<F,TE>>>> TensorReader::readText(
  const std::string &fileName,
  const std::vector<size_t> &lens
) {
  constexpr size_t MAX_BUFFER_SIZE(128*1024*1024);
  std::ifstream stream(fileName.c_str());
  if (stream.fail()) {
    std::stringstream explanation;
    explanation << "Failed to open file \"" << fileName << "\"";
    throw new EXCEPTION(explanation.str());
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
  LOG(1, "TensorReader") << "indexCount=" << A->getElementsCount() << std::endl;
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
    LOG(2, "TensorReader") << "writing " << elementsCount << " values to tensor..." << std::endl;
    A->write(localElementsCount, indices.data(), values.data());
    index += elementsCount;
  }

  return New<AtomicNode<Ptr<Tensor<F,TE>>>>(A);
}

