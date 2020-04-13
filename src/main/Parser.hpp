#ifndef INTERPRETER_DEFINED
#define INTERPRETER_DEFINED

#include <Data.hpp>
#include <algorithms/Algorithm.hpp>
#include <util/LineNumberStream.hpp>
#include <util/SharedPointer.hpp>
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
    Parser(const std::string &fileName);
    ~Parser();

    /**
     * \brief Parses the cc4s algorithms contained in the stream.
     * This method must be called with the same stream content on all processes.
     */
    std::vector<Ptr<Algorithm>> parse();

  protected:
    Ptr<Algorithm> parseAlgorithm();
    std::vector<Argument> parseArguments();
    Argument parseArgument();
    Argument parseImplicitlyNamedArgument();
    Argument parseExplicitlyNamedArgument();
    std::string parseData();
    std::string parseSymbolName();
    Ptr<Data> parseSymbol();
    Ptr<TextData> parseText();
    Ptr<NumericData> parseNumber();
    Ptr<RealData> parseReal(const int64_t sign, const int64_t integerPart);

    void skipIrrelevantCharacters();
    void skipComment();
    void skipWhiteSpaceCharacters();
    void expectCharacter(const char character);

    LineNumberStream stream;
  };
}

#endif

