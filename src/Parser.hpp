#ifndef INTERPRETER_DEFINED
#define INTERPRETER_DEFINED

#include <Data.hpp>
#include <algorithms/Algorithm.hpp>
#include <util/LineNumberStream.hpp>
#include <string>
#include <vector>

namespace cc4s {
  /**
   * \brief Parser for cc4s files specifying the calculation plan, i.e.
   * which algorithms to use in which order.
   */
  class Parser {
  public:
    /**
     * \brief Creates a new interpreter for a cc4s file of the given name.
     * Upon creation the file will be openend but not yet read.
     */
    Parser(std::string const &fileName);
    ~Parser();

    /**
     * \brief Parses the cc4s algorithms contained in the stream.
     * This method must be called with the same stream content on all processes.
     */
    std::vector<Algorithm *> parse();

  protected:
    Algorithm *parseAlgorithm();
    std::vector<Argument> parseArguments();
    Argument parseArgument();
    Argument parseImplicitlyNamedArgument();
    Argument parseExplicitlyNamedArgument();
    std::string parseData();
    std::string parseSymbolName();
    Data *parseSymbol();
    TextData *parseText();
    NumericData *parseNumber();
    RealData *parseReal(int64_t const sign, int64_t const integerPart);

    void skipIrrelevantCharacters();
    void skipComment();
    void skipWhiteSpaceCharacters();
    void expectCharacter(char const character);

    LineNumberStream stream;
  };
}

#endif

