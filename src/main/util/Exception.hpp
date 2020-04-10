/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef EXCEPTION_DEFINED
#define EXCEPTION_DEFINED

#include <string>
#include <sstream>
#include <iostream>

#define EXCEPTION(message) \
  cc4s::DetailedException((message), __FILE__, __LINE__)
#define Assert(condition, message) \
  if (!(condition)) throw new EXCEPTION(message);

namespace cc4s{
  class Exception {
  public:
    virtual std::string getMessage() = 0;
  };

  class DetailedException {
  public:
    DetailedException(
       std::string const &message_, std::string const &file_, int line_
    ): message(message_), file(file_), line(line_), column(0) {
    }
    DetailedException(
       std::string const &message_, std::string const &file_,
        int line_, int column_
    ): message(message_), file(file_), line(line_), column(column_) {
    }
    DetailedException(
       std::stringstream const &stream_, std::string const &file_, int line_
    ): message(stream_.str()), file(file_), line(line_) {
    }
    virtual ~DetailedException() {
    }
    virtual std::string getMessage() {
      std::stringstream sStream;
      sStream << message << std::endl << "\tat " << file << ":" << line;
      if (column > 0) sStream << ":" << column;
      return sStream.str();
    }
  private:
    std::string message, file;
    int line, column;
  };
}

#endif

