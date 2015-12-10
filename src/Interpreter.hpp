#ifndef INTERPRETER_DEFINED
#define INTERPRETER_DEFINED

#include <Data.hpp>
#include <string>

namespace cc4s {
  /**
   * \brief Interpreter for cc4s files specifying the calculation plan, i.e.
   * which algorithms to use in which order.
   */
  class Interpreter {
  public:
    Interpreter();
    ~Interpreter();

    /**
     * \brief Executes the cc4s file with the given name.
     * This method must be called with the same file on all processes.
     */
    void execute(std::string const &fileName);

  protected:
    void parseAlgorithm(std::istream &stream);
    void parseArguments(std::istream &stream);
    std::string parseSymbolName(std::istream &stream);
    TextData *parseText(std::istream &stream);
    NumericData *parseNumber(std::istream &stream);
    RealData *parseReal(
      std::istream &stream, int64_t const sign, int64_t const integerPart
    );
    Data *parseSymbol(std::istream &stream);

    void skipWhiteSpaceCharacters(std::istream &stream);
  };
}

#endif

