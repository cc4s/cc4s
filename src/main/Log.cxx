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

#include <Log.hpp>
#include <Complex.hpp>

#include <sstream>
#include <string>

using namespace cc4s;

int Log::rank(-1);
std::string Log::fileName("cc4s.log");
std::ofstream Log::logStream;
Log::HeaderFunction Log::outHeaderFunction(
  [](const SourceLocation &){ return ""; }
);
Log::HeaderFunction Log::errorHeaderFunction(
  [](const SourceLocation &location){
    std::stringstream header;
    header << location << ": ERROR: ";
    return header.str();
  }
);
Log::HeaderFunction Log::warningHeaderFunction(
  [](const SourceLocation &location){
    std::stringstream header;
    header << location << ": WARNING: ";
    return header.str();
  }
);
Log::HeaderFunction Log::logHeaderFunction(
  [](const SourceLocation &){ return ""; }
);

void Log::setRank(int const rank_) {
  rank = rank_;
  if (rank == 0) {
    logStream.open(
      fileName.c_str(), std::ofstream::out | std::ofstream::trunc
    );
  } else {
    // prevent writing and use this null stream for both: OUT() and LOG()
    logStream.setstate(std::ios_base::badbit);
  }
}

int Log::getRank() {
  return rank;
}

void Log::setFileName(const std::string &fileName_) {
  fileName = fileName_;
}

std::string Log::getFileName() {
  return fileName;
}

