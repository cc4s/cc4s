/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
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
      std::string const &file_, int line_
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
    int getLine() const {
      return line;
    }
  protected:
    std::string file;
    int line;
  };

  inline std::ostream &operator <<(
    std::ostream &stream, SourceLocation const &location
  ) {
    return stream << location.getFile()  << "(" << location.getLine() << ")";
  }
}

#define SOURCE_LOCATION SourceLocation(__FILE__, __LINE__)

#endif
