#include "Log.hpp"

#include <math/Complex.hpp>
#include <sstream>
#include <string>

using namespace cc4s;

int Log::rank(-1);
std::string Log::fileName("cc4s.log");
std::ofstream Log::stream;
Log::HeaderFunction Log::outHeaderFunction(
  [](const SourceLocation &){ return ""; }
);
Log::HeaderFunction Log::errorHeaderFunction(
  [](const SourceLocation &location){
    std::stringstream header;
    header << location << ": \033[1;31mERROR:\033[0m ";
    return header.str();
  }
);
Log::HeaderFunction Log::warningHeaderFunction(
  [](const SourceLocation &location){
    std::stringstream header;
    header << location << ": \033[1;35mWARNING:\033[0m ";
    return header.str();
  }
);
Log::HeaderFunction Log::logHeaderFunction(
  [](const SourceLocation &){ return ""; }
);

void Log::setRank(int const rank_) {
  rank = rank_;
  if (rank == 0) {
    stream.open(
      fileName.c_str(), std::ofstream::out | std::ofstream::trunc
    );
  } else {
    // prevent writing to stdout and to log file
    stream.setstate(std::ios_base::badbit);
    std::cout.setstate(std::ios_base::badbit);
  }
}

int Log::getRank() {
  return rank;
}

void Log::setFileName(const std::string &fileName_) {
  fileName = fileName_;
}

std::string Log::getFileName() {
  return fileName;
}

