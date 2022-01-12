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

#include <algorithms/TensorWriter.hpp>

#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <Real.hpp>
#include <Complex.hpp>
#include <Emitter.hpp>

#include <vector>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(TensorWriter)

/**
 * \brief Writer for tensor data.
 */
Ptr<MapNode> TensorWriter::run(const Ptr<MapNode> &arguments) {
  auto tensor(arguments->getMap("tensor"));
  auto scalarType(tensor->getValue<std::string>("scalarType"));
  auto fileName(arguments->getValue<std::string>("fileName"));
  auto binary(arguments->getValue<bool>("binary", false));

  // if fileName contains '/' change directory
  char currentDirectory[PATH_MAX];
  getcwd(currentDirectory, sizeof(currentDirectory));
  // TODO: support other filesystems
  auto slashPosition(fileName.rfind('/'));
  if (slashPosition != std::string::npos) {
    auto fileDirectory(fileName.substr(0, slashPosition));
    chdir(fileDirectory.c_str());
  }
  auto dotPosition(fileName.rfind('.'));
  ASSERT(
    dotPosition != std::string::npos,
    "'fileName' must contain an extension e.g. '.yaml'"
  );
  auto baseName(fileName.substr(slashPosition+1, dotPosition-slashPosition-1));

  // prepare node containing meta data for writing tensor
  auto writtenTensor(New<MapNode>(SOURCE_LOCATION));
  writtenTensor->setValue("version", std::string("v1.0"));
  writtenTensor->get("scalarType") = tensor->get("scalarType");
  writtenTensor->get("indices") = tensor->get("indices");
  writtenTensor->setValue("binary", binary);

  if (tensor->get("components")) {
    auto writtenComponents(New<MapNode>(SOURCE_LOCATION));
    auto components(tensor->getMap("components"));
    for (auto key: components->getKeys()) {
      auto component(components->getMap(key));
      auto writtenComponent(New<MapNode>(SOURCE_LOCATION));
      writtenComponent->get("data") = component->get("data");
      auto componentBaseName(baseName + "." + key);
      writeData(
        writtenComponent, componentBaseName,
        scalarType, binary, arguments->sourceLocation
      );
      writtenComponents->get(key) = writtenComponent;
    }
    writtenTensor->get("components") = writtenComponents;
  } else {
    writtenTensor->get("dimensions") = tensor->get("dimensions");
    writtenTensor->get("data") = tensor->get("data");
    writeData(
      writtenTensor, baseName, scalarType,
      binary, arguments->sourceLocation
    );
  }

  // change back to current directory
  chdir(currentDirectory);

  // write meta data to yaml file
  Emitter(fileName).emit(writtenTensor);

  // create empty result
  auto result(New<MapNode>(tensor->sourceLocation));
  result->get("tensor") = writtenTensor;
  return result;
}

void TensorWriter::writeData(
  const Ptr<MapNode> &tensorNode,
  const std::string &baseName,
  const std::string &scalarType,
  const bool binary,
  const SourceLocation &sourceLocation
) {
  std::string fileName;
  // multiplex calls to template methods depending on tensor engine and type
  if (binary) {
    fileName = baseName + ".bin";
    if (Cc4s::options->dryRun) {
      using TE = DefaultDryTensorEngine;
      if (scalarType == "real64") {
        writeBinary(
          tensorNode,
          tensorNode->getPtr<TensorExpression<Real<64>,TE>>("data"),
          fileName, sourceLocation
        );
      } else if (scalarType == "complex64") {
        writeBinary(
          tensorNode,
          tensorNode->getPtr<TensorExpression<Complex<64>,TE>>("data"),
          fileName, sourceLocation
        );
      } else {
        ASSERT_LOCATION(
          false, "scalar type '" + scalarType + "' not supported", sourceLocation
        );
      }
    } else {
      using TE = DefaultTensorEngine;
      if (scalarType == "real64") {
        writeBinary(
          tensorNode,
          tensorNode->getPtr<TensorExpression<Real<64>,TE>>("data"),
          fileName, sourceLocation
        );
      } else if (scalarType == "complex64") {
        writeBinary(
          tensorNode,
          tensorNode->getPtr<TensorExpression<Complex<64>,TE>>("data"),
          fileName, sourceLocation
        );
      } else {
        ASSERT_LOCATION(
          false, "scalar type '" + scalarType + "' not supported", sourceLocation
        );
      }
    }
  } else {
    fileName = baseName + ".dat";
    if (Cc4s::options->dryRun) {
      using TE = DefaultDryTensorEngine;
      if (scalarType == "real64") {
        writeText(
          tensorNode,
          tensorNode->getPtr<TensorExpression<Real<64>,TE>>("data"),
          fileName, sourceLocation
        );
      } else if (scalarType == "complex64") {
        writeText(
          tensorNode,
          tensorNode->getPtr<TensorExpression<Complex<64>,TE>>("data"),
          fileName, sourceLocation
        );
      } else {
        ASSERT_LOCATION(
          false, "scalar type '" + scalarType + "' not supported", sourceLocation
        );
      }
    } else {
      using TE = DefaultTensorEngine;
      if (scalarType == "real64") {
        writeText(
          tensorNode,
          tensorNode->getPtr<TensorExpression<Real<64>,TE>>("data"),
          fileName, sourceLocation
        );
      } else if (scalarType == "complex64") {
        writeText(
          tensorNode,
          tensorNode->getPtr<TensorExpression<Complex<64>,TE>>("data"),
          fileName, sourceLocation
        );
      } else {
        ASSERT_LOCATION(
          false, "scalar type '" + scalarType + "' not supported", sourceLocation
        );
      }
    }
  }

  tensorNode->setValue("data", fileName);
}

template <typename F, typename TE>
void TensorWriter::writeText(
  const Ptr<MapNode> &tensorNode,
  const Ptr<TensorExpression<F,TE>> &tensorExpression,
  const std::string &fileName,
  const SourceLocation &sourceLocation
) {
  constexpr size_t MAX_BUFFER_SIZE(128*1024*1024);
  std::ofstream stream(fileName.c_str());
  if (stream.fail()) {
    std::stringstream explanation;
    explanation << "Failed to open file \"" << fileName << "\"";
    throw New<Exception>(explanation.str(), sourceLocation);
  }
  // TODO: use Stream Printer rather than <<
  stream << setprecision(17);

  auto tensor(tensorExpression->evaluate());
  // build dimensions meta data from tensor length if not already present
  if (!tensorNode->get("dimensions")) {
    auto dimensions(New<MapNode>(SOURCE_LOCATION));
    for (size_t d(0); d < tensor->lens.size(); ++d) {
      auto dimension(New<MapNode>(SOURCE_LOCATION));
      dimension->setValue("length", tensor->lens[d]);
      dimension->setValue("length", tensor->dimensions[d]->name);
      dimensions->get(d) = dimension;
    }
    tensorNode->get("dimensions") = dimensions;
  }

  // write the values only on root, all others still pariticipate calling MPI
  const size_t bufferSize(std::min(tensor->getElementsCount(), MAX_BUFFER_SIZE));
  size_t localBufferSize(Cc4s::world->getRank() == 0 ? bufferSize : 0);
  std::vector<size_t> indices(localBufferSize);
  std::vector<F> values(localBufferSize);

  OUT() << "Writing to text file " << fileName << std::endl;
  if (Cc4s::options->dryRun) {
    stream << "# dry-run: no data written" << std::endl;
    return;
  }

  size_t index(0);
  LOG() << "indexCount: " << tensor->getElementsCount() << std::endl;
  while (index < tensor->getElementsCount()) {
    size_t elementsCount(std::min(bufferSize, tensor->getElementsCount()-index));
    size_t localElementsCount(Cc4s::world->getRank() == 0 ? elementsCount : 0);
    for (size_t i(0); i < localElementsCount; ++i) {
      indices[i] = index+i;
    }
    LOG() << "reading " << elementsCount << " values from tensor..." << std::endl;
    tensor->read(localElementsCount, indices.data(), values.data());
    for (size_t i(0); i < localElementsCount; ++i) {
      stream << values[i] << "\n";
    }
    // wait until all processes finished writeing this buffer into the tensor
    Cc4s::world->barrier();
    index += elementsCount;
  }
  LOG() << "Written " << tensor->getElementsCount() <<
    " elements to text file " << fileName << std::endl;
}

template <typename F, typename TE>
void TensorWriter::writeBinary(
  const Ptr<MapNode> &tensorNode,
  const Ptr<TensorExpression<F,TE>> &tensorExpression,
  const std::string &fileName,
  const SourceLocation &sourceLocation
) {
  // open the file
  MPI_File file;
  int mpiError(
    MPI_File_open(
      Cc4s::world->getComm(), fileName.c_str(),
      MPI_MODE_CREATE | MPI_MODE_WRONLY,
      MPI_INFO_NULL, &file
    )
  );
  ASSERT_LOCATION(
    !mpiError, std::string("Failed to open file '") + fileName + "'",
    sourceLocation
  )

  auto tensor(tensorExpression->evaluate());
  if (!tensorNode->get("dimensions")) {
    // build dimensions meta data from tensor length
    auto dimensions(New<MapNode>(SOURCE_LOCATION));
    for (size_t d(0); d < tensor->lens.size(); ++d) {
      auto dimension(New<MapNode>(SOURCE_LOCATION));
      dimension->setValue("length", tensor->lens[d]);
      dimensions->get(d) = dimension;
      // TODO: how to determine dimension type?
    }
    tensorNode->get("dimensions") = dimensions;
  }

  OUT() << "Writing to binary file " << fileName << std::endl;
  if (Cc4s::options->dryRun) return;

  // write tensor data with values from file
  tensor->readToFile(file);
  LOG() << "Written " << sizeof(F)*tensor->getElementsCount() <<
    " bytes to binary file " << fileName << std::endl;

  // done
  MPI_File_close(&file);
}

