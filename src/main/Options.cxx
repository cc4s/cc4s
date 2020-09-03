/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Options.hpp>
#include <string>
#include <sstream>

using namespace cc4s;

Options::Options(int argumentCount, char **arguments) {
  name = "cc4s";
  inFile = "";
  dryRun = false;
  for (int i(0); i < argumentCount; ++i) {
    std::string argument(arguments[i]);
    // deprecated: old options
    if (argument == "-i" || argument == "-file") {
      inFile = arguments[++i];
    } else if (argument == "-n" || argument == "--name") {
      name = arguments[++i];
    } else if (argument == "--dry-run" || argument == "-dryRun") {
      dryRun = true;
    }
  }
  if (inFile == "") inFile = name + ".in";
}

