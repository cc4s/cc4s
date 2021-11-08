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

#ifndef EXCEPTION_DEFINED
#define EXCEPTION_DEFINED

#include <util/SourceLocation.hpp>
#include <util/SharedPointer.hpp>

#include <string>
#include <sstream>
#include <iostream>

#define THROW(message) \
  throw New<cc4s::Exception>((message), SourceLocation(__FILE__, __LINE__))
#define THROW_LOCATION(message, sourceLocation) \
  throw New<cc4s::Exception>((message), (sourceLocation))
#define ASSERT(condition, message) \
  if (!(condition)) THROW((message))
#define ASSERT_LOCATION(condition, message, sourceLocation) \
  if (!(condition)) throw New<Exception>( \
    (message), SourceLocation(__FILE__, __LINE__), \
    New<Exception>("Request from here", (sourceLocation)) \
  );

namespace cc4s{
  class Exception: public std::exception {
  public:
    Exception(
      std::string const &message_, const SourceLocation &sourceLocation_,
      const Ptr<Exception> &cause_ = nullptr
    ): message(message_), sourceLocation(sourceLocation_), cause(cause_) {
    }
    const char *what() const noexcept override {
      return message.c_str();
    }
    SourceLocation getSourceLocation() const noexcept {
      return sourceLocation;
    }
    Ptr<Exception> getCause() const noexcept {
      return cause;
    }
  private:
    std::string message;
    SourceLocation sourceLocation;
    Ptr<Exception> cause;
  };
}

#endif

