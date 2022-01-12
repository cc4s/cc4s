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

#ifndef SOURCE_LOCATION_DEFINED
#define SOURCE_LOCATION_DEFINED

#include <string>
#include <sstream>
#include <ostream>

namespace cc4s {
  class SourceLocation;
  std::ostream &operator <<(std::ostream &stream, SourceLocation const &l);

  class SourceLocation {
  public:
    SourceLocation(): file("<unknown>"), line(0) {
    }
    SourceLocation(SourceLocation const &l): file(l.file), line(l.line) {
    }
    SourceLocation(
      std::string const &file_, size_t line_
    ): file(file_), line(line_) {
    }
    std::string getLocation() const {
      std::stringstream stream;
      stream << *this;
      return stream.str();
    }
    std::string getFile() const {
      return file;
    }
    size_t getLine() const {
      return line;
    }
    bool isValid() const {
      return file != "<unknown>";
    }
  protected:
    std::string file;
    size_t line;
  };

  inline std::ostream &operator <<(
    std::ostream &stream, SourceLocation const &location
  ) {
    return stream << location.getFile()  << "(" << location.getLine() << ")";
  }
}

#define SOURCE_LOCATION SourceLocation(__FILE__, __LINE__)

#endif
