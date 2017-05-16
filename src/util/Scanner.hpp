/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SCANNER_DEFINED
#define SCANNER_DEFINED

#include <math/Complex.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <sstream>
#include <istream>
#include <cstring>

namespace cc4s {
  class Scanner {
  protected:
    static unsigned int constexpr BUFFER_SIZE = 128*1024*1024;
    static unsigned int constexpr REFILL_SIZE = 128;

    std::istream *stream;
    char *buffer, *pos, *end;
    bool eof;

    void refillBuffer() {
      if (pos + REFILL_SIZE > buffer + BUFFER_SIZE && !eof) {
        // move remaining text to beginning of buffer
        std::memmove(buffer, pos, end-pos);
        // set the buffer end to the end of the remaining text
        end = buffer + (end - pos);
        // the next character to read is now at the beginning
        pos = buffer;
        // fill the rest of the buffer from the file
        int64_t size(buffer+BUFFER_SIZE-end);
        int64_t count(stream->read(end, size).gcount());
        LOG(2, "Scanner") << count << " bytes fetched." << std::endl;
        // account the read characters to the buffer
        end += count;
        // if not all requested characters could be read the file is done
        if (count < size) {
          eof = true;
          *end = 0;
        }
      }
    }

  public:
    Scanner(
      std::istream *stream_
    ):
      stream(stream_),
      buffer(new char[BUFFER_SIZE]),
      pos(buffer+BUFFER_SIZE-1), end(pos),
      eof(false)
    {
    }

    std::string nextLine(char const delimiter = '\n') {
      std::stringstream lineStream;
      bool terminated(false);
      do {
        refillBuffer();
        switch (*pos) {
        case 0:
          terminated = true;
          break;
        default:
          char c(*pos++);
          if (c != delimiter) lineStream.put(c);
          else terminated = true;
          break;
        }
      } while (!terminated);
      return lineStream.str();
    }

    template <typename NumberType>
    friend class NumberScanner;
  };

  template <typename NumberType=double>
  class NumberScanner {
  };
  template <>
  class NumberScanner<double> {
    public:
    NumberScanner(Scanner *scanner_): scanner(scanner_) {
    }
    double nextNumber() {
      scanner->refillBuffer();
      return std::strtod(scanner->pos, &scanner->pos);
    }
  protected:
    Scanner *scanner;
  };
  template <>
  class NumberScanner<complex> {
    public:
    NumberScanner(Scanner *scanner_): scanner(scanner_) {
    }
    complex nextNumber() {
      scanner->refillBuffer();
      while (isspace(*scanner->pos) || *scanner->pos == '(') ++scanner->pos;
      // read real part
      double r(std::strtod(scanner->pos, &scanner->pos));
      // skip ','
      ++scanner->pos;
      // read imaginary part
      double i(std::strtod(scanner->pos, &scanner->pos));
      // skip ')'
      ++scanner->pos;
      return complex(r, i);
    }
  protected:
    Scanner *scanner;
  };
}
#endif

