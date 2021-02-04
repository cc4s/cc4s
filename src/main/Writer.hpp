/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef WRITER_DEFINED
#define WRITER_DEFINED

#include <Data.hpp>
#include <Emitter.hpp>
#include <Cc4s.hpp>
#include <tcc/Tensor.hpp>
#include <tcc/TensorRecipe.hpp>
#include <util/SharedPointer.hpp>
#include <string>

namespace cc4s {
  /**
   * \brief Writer for writing entire node states persistently to files,
   * including tensor data.
   * In contrast, the Emitter only writes memory references of tensors,
   * which can only be used during the execution of cc4s.
   */
  class Writer {
  public:
    Writer(
      const std::string &fileName_, const bool useBinary_ = false
    ): fileName(fileName_), useBinary(useBinary_) {
    }
    Ptr<Node> write(const Ptr<Node> &node) {
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
      // base name needed for tensor data
      baseName = fileName.substr(slashPosition+1, dotPosition-slashPosition-1);

      // translate into persistent node tree. Contained tensor data are
      // written as side effects.
      auto persistentNode(persist(node));

      // restore to original directory
      chdir(currentDirectory);

      // finally, write persistent tree to given file
      Emitter(fileName).emit(persistentNode);

      return persistentNode;
    }

  protected:
    // translate into persistent node tree
    Ptr<Node> persist(const Ptr<Node> &node, const std::string &nodePath = "") {
      auto mapNode(node->toMap());
      if (mapNode) return persistMap(mapNode, nodePath);
      auto symbolNode(node->toSymbol());
      if (symbolNode) {
        // symbols are already in persistent form
        return symbolNode;
      } else {
        return persistAtom(node, nodePath);
      }
    }

    Ptr<MapNode> persistMap(
      const Ptr<MapNode> &mapNode,
      const std::string &nodePath
    ) {
      auto persistentMapNode(New<MapNode>(mapNode->sourceLocation));
      for (auto key: mapNode->getKeys()) {
        auto valueNode(mapNode->get(key));
        // NOTE that there may be nullptr keys, if keys have been searched for
        if (valueNode) {
          persistentMapNode->get(key) = persist(valueNode, nodePath + "." + key);
        }
      }
      return persistentMapNode;
    }

    Ptr<Node> persistAtom(
      const Ptr<Node> &node, const std::string &nodePath
    ) {
      // multiplex different tensor types
      Ptr<Node> persistentNode;
      if (!Cc4s::options->dryRun) {
        using TE = DefaultTensorEngine;
        persistentNode = persistTensor<Real<64>,TE>(node, nodePath);
        if (persistentNode) return persistentNode;
        persistentNode = persistTensor<Complex<64>,TE>(node, nodePath);
        if (persistentNode) return persistentNode;
        persistentNode = persistTensor<Real<128>,TE>(node, nodePath);
        if (persistentNode) return persistentNode;
        persistentNode = persistTensor<Complex<128>,TE>(node, nodePath);
        if (persistentNode) return persistentNode;
      } else {
        using TE = DryTensorEngine;
        persistentNode = persistTensor<Real<64>,TE>(node, nodePath);
        if (persistentNode) return persistentNode;
        persistentNode = persistTensor<Complex<64>,TE>(node, nodePath);
        if (persistentNode) return persistentNode;
        persistentNode = persistTensor<Real<128>,TE>(node, nodePath);
        if (persistentNode) return persistentNode;
        persistentNode = persistTensor<Complex<128>,TE>(node, nodePath);
        if (persistentNode) return persistentNode;
      }
      // other types of atomic nodes can be serialized with to_string
      return node;
    }

    template <typename F, typename TE>
    Ptr<MapNode> persistTensor(
      const Ptr<Node> &node,
      const std::string &nodePath
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
      auto persistentTensor(New<MapNode>(SOURCE_LOCATION));
      persistentTensor->setSymbol("type", "tensor");

      // build dimensions from tensor length
      auto dimensions(New<MapNode>(SOURCE_LOCATION));
      for (size_t d(0); d < tensor->lens.size(); ++d) {
        auto dimension(New<MapNode>(SOURCE_LOCATION));
        dimension->setValue<size_t>("length", tensor->lens[d]);
        dimensions->get(d) = dimension;
        // TODO: how to determine dimension type?
      }
      persistentTensor->get("dimensions") = dimensions;
      persistentTensor->setSymbol("scalarType", TypeTraits<F>::getName());

      // write tensor data as side effect
      if (useBinary) {
        writeTensorDataBinary(tensor, nodePath);
      } else {
        writeTensorDataText(tensor, nodePath);
      }

      return persistentTensor;
    }

    template <typename F, typename TE>
    void writeTensorDataText(
      const Ptr<Tensor<F,TE>> &tensor, const std::string &nodePath
    ) {
      constexpr size_t MAX_BUFFER_SIZE(128*1024*1024);
      std::string dataFileName(baseName + nodePath + ".dat");
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
    void writeTensorDataBinary(
      const Ptr<Tensor<F,TE>> &tensor, const std::string &nodePath
    ) {
      std::string dataFileName(baseName + "." + nodePath + ".bin");
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

    std::string fileName, baseName;
    bool useBinary;
  };
}

#endif

