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

#ifndef TCC_OPERATION_DEFINED
#define TCC_OPERATION_DEFINED

#include <tcc/Costs.hpp>
#include <Integer.hpp>

namespace cc4s {
  template <typename TE>
  class Operation {
  public:
    static Natural<128> getFloatingPointOperations() {
      return floatingPointOperations;
    }
    static void addFloatingPointOperations(const Natural<128> ops) {
      floatingPointOperations += ops;
    }

    virtual ~Operation() {
    }

    virtual void execute() = 0;

    virtual size_t getLatestSourceVersion() = 0;

    virtual operator std::string () const = 0;

    /**
     * \brief Costs to evaluate this operation in time and memory
     **/
    Costs costs;

    std::string file;
    size_t line;

  protected:
    /**
     * \brief Total executed floating point operations.
     **/
    static Natural<128> floatingPointOperations;

    Operation(
      const Costs &costs_, const std::string &file_, const size_t line_
    ): costs(costs_), file(file_), line(line_) {
    }
    class ProtectedToken {
    };
  };

  template<typename TE>
  Natural<128> cc4s::Operation<TE>::floatingPointOperations = 0;
}

#endif

