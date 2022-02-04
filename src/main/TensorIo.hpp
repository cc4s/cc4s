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
    static Ptr<TensorDimension> getDimension(const std::string &name);

    /**
     * \brief Static handler routine for writing nodes of the type
     * PointerNode<TensorExpression<F,TE>>. If the given node is of such a type
     * the contained tensor will be written to disk and a MapNode is
     * constructed and returned containing all necessary information to
     * load the written tensor again.
     * nullptr is returned if the given node is not a PointerNode containing
     * a tensor pointer.
     **/
    static Ptr<Node> write(
      const Ptr<Node> &node, const std::string &nodePath, const bool useBinary
    );

    static Ptr<Node> read(
      const Ptr<MapNode> &node, const std::string &nodePath
    );

  protected:
    /**
     * \brief Serialization version. Only objects written by
     * a matching serialization version can be read. Increase this version
     * if the data written is no longer backward compatible.
     **/
    static const Natural<32> VERSION;

    template <typename F, typename TE>
    static Ptr<MapNode> writeTensor(
      const Ptr<Node> &node,
      const std::string &nodePath,
      const bool useBinary
    );

    template <typename F, typename TE>
    static Ptr<PointerNode<Object>> readTensor(
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
    static Ptr<Tensor<F,TE>> readTensorElementsText(
      const std::string &fileName,
      const std::vector<size_t> &lens,
      const Ptr<TensorNonZeroConditions> &nonZeroConditions,
      const SourceLocation &sourceLocation
    );

    template <typename F, typename TE>
    static Ptr<Tensor<F,TE>> readTensorElementsBinary(
      const std::string &fileName,
      const std::vector<size_t> &lens,
      const Ptr<TensorNonZeroConditions> &nonZeroConditions,
      const SourceLocation &sourceLocation
    );

    template <typename F, typename TE>
    static Ptr<Tensor<F,TE>> readTensorElementsBinaryDense(
      const std::string &fileName,
      const std::vector<size_t> &lens,
      const Ptr<TensorNonZeroConditions> &nonZeroConditions,
      const SourceLocation &sourceLocation
    );

    static bool WRITE_REGISTERED, READ_REGISTERED;
  };
}

#endif

