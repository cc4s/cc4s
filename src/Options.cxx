/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Options.hpp>
#include <string>
#include <sstream>

using namespace cc4s;

Options::Options(int argumentCount, char **arguments) {
  profile = DEFAULT_PROFILE;
  storeV = DEFAULT_STORE_V;
  stridedIo = DEFAULT_STRIDED_IO;
  file = "calculation.cc4s";
  logLevel = DEFAULT_LOG_LEVEL;
  logFile = "cc4s.log";
  rank = DEFAULT_RANK;
  accuracy = DEFAULT_ACCURACY;
  dryRun = false;
  for (int i(0); i < argumentCount; ++i) {
    std::string argument(arguments[i]);
    if (argument == "-file" || argument == "-i") {
      file = arguments[++i];
    } else if (argument == "-profile") { // TODO: remove?
      std::stringstream stream(arguments[++i]);
      stream >> profile;
    } else if (argument == "-storeV") { // TODO: remove?
      std::stringstream stream(arguments[++i]);
      stream >> storeV;
    } else if (argument == "-stridedIo") { // TODO: remove?
      std::stringstream stream(arguments[++i]);
      stream >> stridedIo;
    } else if (argument == "-logLevel") {
      std::stringstream stream(arguments[++i]);
      stream >> logLevel;
    } else if (argument == "--yaml") {
      yamlFile = arguments[++i];
    } else if (argument == "-logFile" || argument == "-o") {
      logFile = arguments[++i];
    } else if (argument == "-rank") { // TODO: remove?
      std::stringstream stream(arguments[++i]);
      stream >> rank;
    } else if (argument == "-accuracy") { // TODO: remove?
      std::stringstream stream(arguments[++i]);
      stream >> accuracy;
    } else if (argument == "-dryRun") {
      dryRun = true;
    }
  }
}

