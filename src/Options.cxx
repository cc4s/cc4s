/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include "Options.hpp"
#include <string>
#include <sstream>

Options::Options(int argumentCount, char **arguments) {
  for (int i(0); i < argumentCount; ++i) {
    std::string argument(arguments[i]);
    nG = DEFAULT_NG;
    nv = DEFAULT_NV;
    no = DEFAULT_NO;
    niter = DEFAULT_NITER;
    profile = DEFAULT_PROFILE;

    if (argument == "-nG") {
      std::stringstream stream(arguments[++i]);
      stream >> nG;
      if (nG < 0) nG = DEFAULT_NG;
    } else if (argument == "-no") {
      std::stringstream stream(arguments[++i]);
      stream >> no;
      if (nG < 0) no = DEFAULT_NO;
    } else if (argument == "-nv") {
      std::stringstream stream(arguments[++i]);
      stream >> nv;
      if (nG < 0) nv = DEFAULT_NV;
    } else if (argument == "-niter") {
      std::stringstream stream(arguments[++i]);
      stream >> niter;
      if (nG < 0) niter = DEFAULT_NITER;
    } else if (argument == "-profile") {
      std::stringstream stream(arguments[++i]);
      stream >> profile;
      if (nG < 0) profile = DEFAULT_PROFILE;
    }
  }
}

