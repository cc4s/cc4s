#ifndef TENSOR_IO_DEFINED
#define TENSOR_IO_DEFINED

#include <Writer.hpp>
#include <util/Scanner.hpp>
#include <algorithms/Algorithm.hpp>
#include <tcc/Tensor.hpp>
#include <tcc/TensorRecipe.hpp>

namespace cc4s {
  class TensorIo {
  public:
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
      // not my type...
      return nullptr;
    }

    template <typename F, typename TE>
    static Ptr<MapNode> writeTensor(
      const Ptr<Node> &node,
      const std::string &nodePath,
      const bool useBinary
    ) {
      Ptr<Tensor<F,TE>> tensor;
      // check if node contains tensor or tensor recipe
      auto tensorNode(node->toAtom<Ptr<Tensor<F,TE>>>());
      if (tensorNode) {
        tensor = tensorNode->value;
      } else {
        auto tensorRecipeNode(node->toAtom<Ptr<TensorRecipe<F,TE>>>());
        if (!tensorRecipeNode) return nullptr;
        auto tensorRecipe(tensorRecipeNode->value);
        tensorRecipe->execute();
        tensor = tensorRecipe->getResult();
      }
      auto writtenTensor(New<MapNode>(SOURCE_LOCATION));

      // build dimensions from tensor length
      auto dimensions(New<MapNode>(SOURCE_LOCATION));
      for (size_t d(0); d < tensor->lens.size(); ++d) {
        auto dimension(New<MapNode>(SOURCE_LOCATION));
        dimension->setValue<size_t>("length", tensor->lens[d]);
        dimensions->get(d) = dimension;
        // TODO: how to determine dimension type?
      }
      writtenTensor->get("dimensions") = dimensions;
      writtenTensor->setSymbol("scalarType", TypeTraits<F>::getName());

      // write tensor data as side effect
      if (useBinary) {
        writeTensorDataBinary(tensor, nodePath);
      } else {
        writeTensorDataText(tensor, nodePath);
      }

      return writtenTensor;
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

    static int WRITE_REGISTERED;
  };
}

#endif

