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

#include <Options.hpp>
#include <string>
#include <sstream>

using namespace cc4s;

Options::Options(int argumentCount, char **arguments) {
  name = "cc4s";
  inFile = "";
  dryRun = false;
  for (int i(0); i < argumentCount; ++i) {
    std::string argument(arguments[i]);
    // deprecated: old options
    if (argument == "-i" || argument == "-file") {
      inFile = arguments[++i];
    } else if (argument == "-n" || argument == "--name") {
      name = arguments[++i];
    } else if (argument == "--dry-run" || argument == "-dryRun") {
      dryRun = true;
    }
  }
  if (inFile == "") inFile = name + ".in";
}

