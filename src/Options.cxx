/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "Options.hpp"
#include <string>
#include <sstream>

Options::Options(int argumentCount, char **arguments) {
  nw = DEFAULT_NW;
  niter = DEFAULT_NITER;
  profile = DEFAULT_PROFILE;
  storeV = DEFAULT_STORE_V;
  for (int i(0); i < argumentCount; ++i) {
    std::string argument(arguments[i]);
    if (argument == "-nw") {
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
    }
  }
}

