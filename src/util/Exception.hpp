/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef EXCEPTION_DEFINED
#define EXCEPTION_DEFINED

#include <string>
#include <iostream>

#define Exception(message) DetailedException(message, __FILE__, __LINE__)
#define Assert(condition, message) \
  if (!(condition)) throw new Exception(message);

namespace cc4s{
  class Exception {
    public:
      virtual std::string getMessage() = 0;
  };

  class DetailedException {
    public:
      DetailedException(
         std::string const &message_, std::string const &file_, int line_
      ): message(message_), file(file_), line(line_) {
      }
      virtual std::string getMessage() {
        std::stringstream sstream;
        sstream << message << std::endl << "\tat " << file << " (" << line << ")";
        return sstream.str();
      }
    private:
      std::string message, file;
      int line;
  };
}

#endif

