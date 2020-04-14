/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Options.hpp>
#include <string>
#include <sstream>

using namespace cc4s;

Options::Options(int argumentCount, char **arguments) {
  inFile = "input.yaml";
  outFile = "output.yaml";
  logFile = "cc4s.log";
  logLevel = DEFAULT_LOG_LEVEL;
  dryRun = false;
  for (int i(0); i < argumentCount; ++i) {
    std::string argument(arguments[i]);
    // deprecated: old options
    if (argument == "-i" || argument == "-file") {
      inFile = arguments[++i];
    } else if (argument == "-o" || argument == "-yamlFile") {
      outFile = arguments[++i];
    } else if (argument == "--log-file" || argument == "-logFile") {
      logFile = arguments[++i];
    } else if (argument == "--log-level" || argument == "-logLevel") {
      std::stringstream stream(arguments[++i]);
      stream >> logLevel;
    } else if (argument == "--dry-run" || argument == "-dryRun") {
      dryRun = true;
    }
  }
}

