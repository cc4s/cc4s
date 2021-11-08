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

#include <algorithms/Read.hpp>

#include <Cc4s.hpp>
#include <Reader.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(Read)

/**
 * \brief Interface to Reader.
 */
Ptr<MapNode> Read::run(const Ptr<MapNode> &arguments) {
  auto fileName(arguments->getValue<std::string>("fileName"));

  try {
    auto data(Reader(fileName).read());
    // create result
    auto result(New<MapNode>(data->sourceLocation));
    result->get("data") = data;
    return result;
  } catch (const Ptr<Exception> &cause) {
    throw New<Exception>(
      "Failed to read file '" + fileName + "'",
      arguments->get("fileName")->sourceLocation,
      cause
    );
  }
}

