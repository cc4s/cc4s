/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef SCANNER_DEFINED
#define SCANNER_DEFINED

#include <util/LineNumberStream.hpp>
#include <sstream>
#include <ctf.hpp>

namespace cc4s {
  class Scanner {
  protected:
    static int constexpr LINE_BUFFER_SIZE = 256;
    static int constexpr NUMBER_BUFFER_SIZE = 256;

    LineNumberStream *stream;

    // TODO: factor out common routines with Parser
    void skipWhiteSpaceCharacters() {
      while (isspace(stream->peek())) stream->get();
    }
    void unexpectedCharacter(std::string const &expectedCharacter) {
      std::stringstream sStream;
      sStream <<
        "Expected '" << expectedCharacter << "', got '" << stream->peek() << "'";
      throw new DetailedException(
        sStream.str(), stream->getSource(), stream->getLine(), stream->getColumn()
      );
    }

  public:
    Scanner(
      LineNumberStream *stream_
    ): stream(stream_) {
    }

    std::stringstream nextLine(char const delimiter = '\n') {
      std::stringstream lineStream;
      bool terminated(false);
      do {
        switch (stream->peek()) {
        case EOF:
          terminated = true;
          break;
        default:
          char c(stream->get());
          if (c != delimiter) lineStream.put(c);
          else terminated = true;
          break;
        }
      } while (!terminated);
      return lineStream;
    }

    double nextReal() {
      char buffer[NUMBER_BUFFER_SIZE];
      skipWhiteSpaceCharacters();
      int i(0);
      bool terminated(false);
      do {
        // TODO: cleanup
        switch (stream->peek()) {
        case '0':
        case '1':
        case '2':
        case '3':
        case '4':
        case '5':
        case '6':
        case '7':
        case '8':
        case '9':
        case '+':
        case '-':
        case '.':
        case 'e':
        case 'E':
          buffer[i++] = stream->get();
          break;
        case ' ':
        case '\n':
        case '\t':
          terminated = true;
          break;
        default:
          unexpectedCharacter("digit");
          break;
        }
      } while (!terminated);
      buffer[i] = 0;
      char *end;
      double result(std::strtod(buffer, &end));
      if (*end == 0) {
        return result;
      } else {
        throw new DetailedException(
          "real was expected", stream->getSource(), stream->getLine(), stream->getColumn()
        );
      }
    }
  };
}
#endif

