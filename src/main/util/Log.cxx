/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include "Log.hpp"

#include <math/Complex.hpp>
#include <sstream>
#include <string>

using namespace cc4s;

/*
std::ostream &LogStream::prepare(
  std::string const &sourceFileName,
  const size_t sourceFileLine,
  const unsigned int level,
  const std::string &category
) {
  Time time(Time::getCurrentRealTime());
  time -= startTime;
  logFile << time << " ";
  std::ostream *log(&logFile);
  if (logLevel >= level) {
    for (unsigned int i(0); i < level; ++i) {
      std::cout << indent.c_str();
    }
    std::stringstream fraction;
    fraction << std::fixed << 1.0*(time.getFractions()) / Time::FRACTIONS;
    std::cout << time.getSeconds() << fraction.str().substr(1,4) << " ";
    // next puts should go to logFile and std::out, done by this->put
    log = this;
  }
  if (category.length() > 0) {
    (*log) << category << ": ";
  } else if (sourceFileName.compare(0, 9, "src/main/") == 0) {
    // remove redundant part of source name
    (*log) << sourceFileName.substr(9) << ':' << sourceFileLine << ": ";
  } else {
    // name is not a c file, take full name
    (*log) << sourceFileName << ':' << sourceFileLine << ": ";
  }
  return *log;
}
*/

int Log::rank(-1);
std::string Log::fileName("cc4s.log");
std::ofstream Log::stream;
Log::HeaderFunction Log::outHeaderFunction(
  [](const SourceLocation &){ return ""; }
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

