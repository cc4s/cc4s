#ifndef INTERPRETER_DEFINED
#define INTERPRETER_DEFINED

#include <Data.hpp>
#include <Algorithm.hpp>
#include <string>
#include <vector>

namespace cc4s {
  /**
   * \brief Wrapper for an input stream providing tracking of the current
   * line and the current column during stream reading.
   */
  class LineNumberStream {
  public:
    LineNumberStream(
      std::istream *stream_, int const tabWidth_=2
    ): line(1), column(1), tabWidth(tabWidth_), stream(stream_) {
    }
    ~LineNumberStream() {
      delete stream;
    }
    int peek() {
      return stream->peek();
    }
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
    std::istream *getStream() { return stream; }
    /**
     * \breif Returns the line of the next character to be read.
     */
    int getLine() { return line; }
    /**
     * \breif Returns the column of the next character to be read.
     */
    int getColumn() { return column; }
  protected:
    int line, column, tabWidth;
    std::istream *stream;
  };

  /**
   * \brief Interpreter for cc4s files specifying the calculation plan, i.e.
   * which algorithms to use in which order.
   */
  class Interpreter {
  public:
    /**
     * \brief Creates a new interpreter for a cc4s file of the given name.
     * Upon creation the file will be openend but not yet read.
     */
    Interpreter(std::string const &fileName);
    ~Interpreter();

    /**
     * \brief Executes the cc4s commands.
     * This method must be called with the same stream content on all processes.
     */
    void execute();

  protected:
    void parseAlgorithm();
    std::vector<Argument> parseArguments();
    Argument parseArgument();
    Argument parseImplicitlyNamedArgument();
    Argument parseExplicitlyNamedArgument();
    std::string parseData();
    std::string parseSymbolName();
    TextData *parseText();
    NumericData *parseNumber();
    RealData *parseReal(int64_t const sign, int64_t const integerPart);
    Data *parseSymbol();

    void skipWhiteSpaceCharacters();
    void expectCharacter(char const character);

    LineNumberStream stream;
  };
}

#endif

