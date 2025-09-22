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

#ifndef WRITER_DEFINED
#define WRITER_DEFINED

#include <Node.hpp>
#include <Emitter.hpp>
#include <Cc4s.hpp>
#include <SharedPointer.hpp>
#include <Log.hpp>

#include <string>

namespace cc4s {

  struct FlagsIO {
    bool useBinary;
    bool rowMajor;
  };

  /**
   * \brief Writer for writing entire node states persistently to files,
   * including tensor data.
   * In contrast, the Emitter only writes memory references of tensors,
   * which can only be used during the execution of cc4s.
   */
  class Writer {
  public:
    Writer(
      const std::string &fileName_, const FlagsIO flagsIO_
    ): fileName(fileName_), flagsIO(flagsIO_) {
    }
    Ptr<Node> write(const Ptr<Node> &node);

    using WriteFunction = std::function<
      Ptr<Node>(const Ptr<Node> &node,
                const std::string &nodePath,
                const FlagsIO flagsIO)
      >;
    static int registerWriteFunction(
      const std::string &name, WriteFunction writeFunction
    );

  protected:
    static std::map<std::string, WriteFunction> writeFunctions;

    // translate into persistent node tree
    Ptr<Node> write(const Ptr<Node> &node, const std::string &nodePath);

    Ptr<MapNode> writeMap(
      const Ptr<MapNode> &mapNode,
      const std::string &nodePath
    );

    Ptr<Node> writeAtom(
      const Ptr<Node> &node, const std::string &nodePath
    );

    std::string fileName;
    FlagsIO flagsIO;
  };

}

#endif

