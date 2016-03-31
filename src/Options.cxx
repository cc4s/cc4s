/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Options.hpp>
#include <string>
#include <sstream>

using namespace cc4s;

Options::Options(int argumentCount, char **arguments) {
  nw = DEFAULT_NW;
  niter = DEFAULT_NITER;
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
    if (argument == "-file") {
      file = arguments[++i];
    } else if (argument == "-nw") {
      std::stringstream stream(arguments[++i]);
      stream >> nw;
      if (nw < 0) nw = DEFAULT_NW;
    } else if (argument == "-niter") {
      std::stringstream stream(arguments[++i]);
      stream >> niter;
      if (niter < 0) niter = DEFAULT_NITER;
    } else if (argument == "-profile") {
      std::stringstream stream(arguments[++i]);
      stream >> profile;
    } else if (argument == "-storeV") {
      std::stringstream stream(arguments[++i]);
      stream >> storeV;
    } else if (argument == "-stridedIo") {
      std::stringstream stream(arguments[++i]);
      stream >> stridedIo;
    } else if (argument == "-logLevel") {
      std::stringstream stream(arguments[++i]);
      stream >> logLevel;
    } else if (argument == "-logFile") {
      logFile = arguments[++i];
    } else if (argument == "-rank") {
      std::stringstream stream(arguments[++i]);
      stream >> rank;
    } else if (argument == "-accuracy") {
      std::stringstream stream(arguments[++i]);
      stream >> accuracy;
    } else if (argument == "-dryRun") {
      dryRun = true;
    }
  }
}

