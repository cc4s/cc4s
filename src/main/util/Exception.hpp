/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
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

