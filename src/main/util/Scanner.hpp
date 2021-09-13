#ifndef SCANNER_DEFINED
#define SCANNER_DEFINED

#include <math/Real.hpp>
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
#ifdef DEBUG
        LOG() << count << " bytes fetched." << std::endl;
#endif
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
      // FIXME: replace with vector<char>
      buffer(new char[BUFFER_SIZE]),
      pos(buffer+BUFFER_SIZE-1), end(pos),
      eof(false)
    {
    }

    ~Scanner() {
      delete[] buffer;
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

  // double precision float
  template <typename NumberType=Real<64>>
  class NumberScanner {
  };
  template <>
  class NumberScanner<Real<64>> {
    public:
    NumberScanner(Scanner *scanner_): scanner(scanner_) {
    }
    Real<64> nextNumber() {
      scanner->refillBuffer();
      return scanReal(&scanner->pos);
    }
    static Real<64> scanReal(char **position) {
      return std::strtod(*position, position);
    }
  protected:
    Scanner *scanner;
  };

  // quadruple precision float
  template <>
  class NumberScanner<Real<128>> {
    public:
    NumberScanner(Scanner *scanner_): scanner(scanner_) {
    }
    Real<128> nextNumber() {
      scanner->refillBuffer();
      return scanReal(&scanner->pos);
    }
    static Real<128> scanReal(char **position) {
      return strtoflt128(*position, position);
    }
  protected:
    Scanner *scanner;
  };

  // complex 64 bit
  template <>
  class NumberScanner<Complex<64>> {
    public:
    NumberScanner(Scanner *scanner_): scanner(scanner_) {
    }
    Complex<64> nextNumber() {
      scanner->refillBuffer();
      while (isspace(*scanner->pos) || *scanner->pos == '(') ++scanner->pos;
      // read real part
      Real<64> r(NumberScanner<Real<64>>::scanReal(&scanner->pos));
      // skip ','
      ++scanner->pos;
      // read imaginary part
      Real<64> i(NumberScanner<Real<64>>::scanReal(&scanner->pos));
      // skip ')'
      ++scanner->pos;
      return Complex<64>(r, i);
    }
  protected:
    Scanner *scanner;
  };

  // complex 128 bit
  template <>
  class NumberScanner<Complex<128>> {
    public:
    NumberScanner(Scanner *scanner_): scanner(scanner_) {
    }
    Complex<128> nextNumber() {
      scanner->refillBuffer();
      while (isspace(*scanner->pos) || *scanner->pos == '(') ++scanner->pos;
      // read real part
      Real<128> r(NumberScanner<Real<128>>::scanReal(&scanner->pos));
      // skip ','
      ++scanner->pos;
      // read imaginary part
      Real<128> i(NumberScanner<Real<128>>::scanReal(&scanner->pos));
      // skip ')'
      ++scanner->pos;
      return Complex<128>(r, i);
    }
  protected:
    Scanner *scanner;
  };
}
#endif

