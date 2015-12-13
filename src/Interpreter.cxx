#include <Interpreter.hpp>

#include <Algorithm.hpp>
#include <util/Exception.hpp>
#include <fstream>
#include <locale>

using namespace cc4s;

// TODO: create ParseException and LineNumberStream

Interpreter::Interpreter(
  std::string const &fileName
): stream(new std::ifstream(fileName)) {
  std::ifstream *fileStream(dynamic_cast<std::ifstream *>(stream.getStream()));
  if (!fileStream->is_open()) {
    std::stringstream sStream;
    sStream << "Failed to open file " << fileName;
    throw new Exception(sStream.str());
  }
}

Interpreter::~Interpreter() {
}

void Interpreter::execute() {
  skipWhiteSpaceCharacters();
  while (stream.peek() > 0) {
    parseAlgorithm();
    skipWhiteSpaceCharacters();
  }
}

void Interpreter::parseAlgorithm() {
  // an algorithm starts with the name
  std::string algorithmName(parseSymbolName());
  std::function<Algorithm *(std::vector<Argument>)> createFunction(
    Algorithm::get(algorithmName)
  );
  skipWhiteSpaceCharacters();
  // next comes its input arguments
  std::vector<Argument> arguments(parseArguments());
  skipWhiteSpaceCharacters();
  // next comes its output arguments
  std::vector<Argument> outputArguments(parseArguments());
  skipWhiteSpaceCharacters();
  // the algorithm must end with a period
  expectCharacter('.');
  // currently there is no distinction between input and output arguments
  arguments.insert(
    arguments.end(), outputArguments.begin(), outputArguments.end()
  );
  // create an instance of the algorithm
  Algorithm *algorithm(createFunction(arguments));
  // execute the algorithm
  algorithm->run();
  delete algorithm;
}

std::vector<Argument> Interpreter::parseArguments() {
  std::vector<Argument> arguments;
  expectCharacter('[');
  skipWhiteSpaceCharacters();
  while (stream.peek() != ']') {
    arguments.push_back(parseArgument());
    skipWhiteSpaceCharacters();
  }
  expectCharacter(']');
  return arguments;
}

Argument Interpreter::parseArgument() {
  if (stream.peek() == '(') return parseExplicitlyNamedArgument();
  else return parseImplicitlyNamedArgument();
}

Argument Interpreter::parseImplicitlyNamedArgument() {
  std::string argumentName(parseSymbolName());
  return Argument(argumentName, argumentName);
}

Argument Interpreter::parseExplicitlyNamedArgument() {
  // first character must be '('
  stream.get();
  skipWhiteSpaceCharacters();
  std::string argumentName(parseSymbolName());
  skipWhiteSpaceCharacters();
  std::string dataName(parseSymbolName());
  skipWhiteSpaceCharacters();
  expectCharacter(')');
  return Argument(argumentName, dataName);
}

std::string Interpreter::parseData() {
  char character(stream.peek());
  if (isalpha(character)) {
    return parseSymbolName();
  } else if (isdigit(character)) {
    return parseNumber()->getName();
  } else if (character == '"') {
    return parseText()->getName();
  } else {
    throw new Exception("Constant or symbol expression expected");
  }
}

std::string Interpreter::parseSymbolName() {
  std::stringstream sStream;
  // the first character is expected to be a alphabetic character
  sStream.put(stream.get());
  char c;
  while (isalpha(c = stream.peek()) || isdigit(c)) {
    sStream.put(stream.get());
  }
  return sStream.str();
}

TextData *Interpreter::parseText() {
  std::stringstream sStream;
  // TODO: parse escape sequences
  // the first character is expected to be a double quote '"'
  stream.get();
  char c;
  while ((c = stream.peek()) > 0 && c != '"') {
    sStream.put(stream.get());
  }
  expectCharacter('"');
  return new TextData(sStream.str());
}

NumericData *Interpreter::parseNumber() {
  // the first character can be a sign
  int64_t sign(1);
  switch (stream.peek()) {
  case '-':
    sign = -1;
  case '+':
    stream.get();
  }
  // the next character must be a digit
  if (!isdigit(stream.peek())) throw new Exception("Digit expected");
  int64_t integer(stream.get() - '0');
  while (isdigit(stream.peek())) {
    integer *= 10;
    integer += stream.get() - '0';
  }
  if (stream.peek() == '.') return parseReal(sign, integer);
  else return new IntegerData(sign * integer);
}

RealData *Interpreter::parseReal(int64_t const sign, int64_t const integerPart){
  // the first character is expected to be the decimal point
  stream.get();
  int64_t numerator(0), denominator(1);
  while (isdigit(stream.peek())) {
    numerator *= 10; denominator *= 10;
    numerator += stream.get() - '0';
  }
  return new RealData(sign * (integerPart + double(numerator) / denominator));
}

Data *Interpreter::parseSymbol() {
  std::string symbolName(parseSymbolName());
  Data *symbol(Data::get(symbolName));
  return symbol ? symbol : new Data(symbolName);
}


void Interpreter::skipWhiteSpaceCharacters() {
  while (isspace(stream.peek())) stream.get();
}

void Interpreter::expectCharacter(char const character) {
  if (stream.get() != character) {
    std::stringstream sStream;
    sStream << "Expected '" << character << "'";
    throw new Exception(sStream.str());
  }
}

