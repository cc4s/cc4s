#include <Interpreter.hpp>

#include <fstream>
#include <locale>

using namespace cc4s;

// TODO: create ParseException and LineNumberStream

Interpreter::Interpreter() {
}

Interpreter::~Interpreter() {
}

void Interpreter::execute(std::string const &fileName) {
  std::ifstream stream(fileName);
  skipWhiteSpaceCharacters(stream);
  while (stream.peek() > 0) {
    // an algorithm starts with the name
  }
}

void Interpreter::parseAlgorithm(std::istream &stream) {
  
}

void Interpreter::parseArguments(std::istream &stream) {
}

std::string Interpreter::parseSymbolName(std::istream &stream) {
  std::stringstream sStream;
  // the first character is expected to be a alphabetic character
  sStream.put(stream.get());
  char c;
  while (isalpha(c = stream.peek()) || isdigit(c)) {
    sStream.put(stream.get());
  }
  return sStream.str();
}

TextData *Interpreter::parseText(std::istream &stream) {
  std::stringstream sStream;
  // TODO: parse escape sequences
  // the first character is expected to be a double quote '"'
  stream.get();
  char c;
  while ((c = stream.peek()) > 0 && c != '"') {
    sStream.put(stream.get());
  }
  if (c != '"') throw new Exception("ending \" expected");
  stream.get();
  return new TextData(sStream.str());
}

NumericData *Interpreter::parseNumber(std::istream &stream) {
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
  if (stream.peek() == '.') return parseReal(stream, sign, integer);
  else return new IntegerData(sign * integer);
}

RealData *Interpreter::parseReal(
  std::istream &stream, int64_t const sign, int64_t const integerPart
) {
  // the first character is expected to be the decimal point
  stream.get();
  int64_t numerator(0), denominator(1);
  while (isdigit(stream.peek())) {
    numerator *= 10; denominator *= 10;
    numerator += stream.get() - '0';
  }
  return new RealData(sign * (integerPart + double(numerator) / denominator));
}

Data *Interpreter::parseSymbol(std::istream &stream) {
  std::string symbolName(parseSymbolName(stream));
  Data *symbol(Data::get(symbolName));
  return symbol ? symbol : new Data(symbolName);
}


void Interpreter::skipWhiteSpaceCharacters(std::istream &stream) {
  while (isspace(stream.peek())) stream.get();
}

