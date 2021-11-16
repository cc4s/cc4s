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

#ifndef TENSOR_WRITER_DEFINED 
#define TENSOR_WRITER_DEFINED

#include <algorithms/Algorithm.hpp>
#include <tcc/TensorExpression.hpp>
#include <util/SourceLocation.hpp>

namespace cc4s {
  class TensorWriter: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(TensorWriter)
    Ptr<MapNode> run(const Ptr<MapNode> &arguments) override;
  protected:
    void writeData(
      const Ptr<MapNode> &tensorNode,
      const std::string &fileName,
      const std::string &scalarType,
      const bool binary,
      const SourceLocation &sourceLocation
    );
    template <typename F, typename TE>
    void writeText(
      const Ptr<MapNode> &tensorNode,
      const Ptr<TensorExpression<F,TE>> &tensorExpression,
      const std::string &fileName,
      const SourceLocation &sourceLocation
    );
    template <typename F, typename TE>
    void writeBinary(
      const Ptr<MapNode> &tensorNode,
      const Ptr<TensorExpression<F,TE>> &tensorExpression,
      const std::string &fileName,
      const SourceLocation &sourceLocation
    );
  };
}

#endif

