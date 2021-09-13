#ifndef OPTIONS_DEFINED
#define OPTIONS_DEFINED

#include <string>

namespace cc4s {
  struct Options {

    std::string inFile, name;
    int logLevel;
    bool dryRun;

    const int DEFAULT_LOG_LEVEL = 1;

    Options(int argumentCount, char **arguments);
  };
}

#endif

