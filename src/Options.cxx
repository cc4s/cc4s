/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Options.hpp>
#include <string>
#include <sstream>

using namespace cc4s;

Options::Options(int argumentCount, char **arguments) {
  file = "calculation.cc4s";
  logLevel = DEFAULT_LOG_LEVEL;
  logFile = "cc4s.log";
  yamlFile = "cc4s.yaml";
  dryRun = false;
  for (int i(0); i < argumentCount; ++i) {
    std::string argument(arguments[i]);
    // deprecated: old options
    if (argument == "-i" || argument == "-file") {
      file = arguments[++i];
    } else if (argument == "--log-level" || argument == "-logLevel") {
      std::stringstream stream(arguments[++i]);
      stream >> logLevel;
    } else if (argument == "--yaml") {
      yamlFile = arguments[++i];
    } else if (argument == "-o" || argument == "-logFile") {
      logFile = arguments[++i];
    } else if (argument == "--dry-run" || argument == "-dryRun") {
      dryRun = true;
    }
  }
}

