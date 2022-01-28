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

#include <algorithms/Write.hpp>

#include <Cc4s.hpp>
#include <Writer.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(Write)

/**
 * \brief Interface to Writer.
 */
Ptr<MapNode> Write::run(const Ptr<MapNode> &arguments) {
  auto source(arguments->get("source"));
  ASSERT_LOCATION(source, "expecting key 'source'", arguments->sourceLocation);
  auto fileName(arguments->getValue<std::string>("fileName"));
  // TODO: use options
  auto useBinary(arguments->getValue<bool>("binary", false));

  auto persistentSource(Writer(fileName, useBinary).write(source));

  // create result
  auto result(New<MapNode>(source->sourceLocation));
  result->get("persitentSource") = persistentSource;
  return result;
}

