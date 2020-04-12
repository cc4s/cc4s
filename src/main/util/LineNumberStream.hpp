#ifndef LINE_NUMBER_STREAM_DEFINED
#define LINE_NUMBER_STREAM_DEFINED

#include <util/SharedPointer.hpp>
#include <string>

namespace cc4s {
  /**
   * \brief Wrapper for an input stream providing tracking of the current
   * line and the current column during stream reading.
   */
  class LineNumberStream {
  public:
    /**
     * \brief Creates a wrapper for the given stream.
     * The wrapper takes ownership of the given pointer.
     */
    LineNumberStream(
      const PTR(std::istream) &stream_,
      const std::string &source_, const size_t tabWidth_=2
    ):
      stream(stream_), line(1), column(1), tabWidth(tabWidth_), source(source_)
    { }
    /**
     * \brief Peeks one character from the underlying stream.
     */
    int peek() {
      return stream->peek();
    }
    /**
     * \brief Reads and returns one character from the underlying stream
     * while keeping track of the line and column number.
     */
    int get() {
      int character(stream->get());
      switch (character) {
      case '\n':
        ++line; column = 1; break;
      case '\t':
        column += tabWidth; break;
      default:
        ++column;
      }
      return character;
    }
    /**
     * \brief Returns the underlying std::istream.
     */
    PTR(std::istream) getStream() { return stream; }
    /**
     * \brief Returns the source name of the underlying stream,
     * usually its file name.
     */
    std::string getSource() { return source; }
    /**
     * \brief Returns the line of the next character to be read.
     */
    int getLine() { return line; }
    /**
     * \brief Returns the column of the next character to be read.
     */
    int getColumn() { return column; }

  protected:
    PTR(std::istream) stream;
    size_t line, column, tabWidth;
    std::string source;
  };
}

#endif

