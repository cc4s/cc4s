/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include "Log.hpp"

#include <math/Complex.hpp>
#include <sstream>
#include <string>

using namespace cc4s;

LogStream::LogStream(
  std::string const &logFileName,
  const unsigned int logLevel_,
  std::string const &indent_
):
  std::ostream(&logBuffer),
  logFile(logFileName.c_str(), std::ofstream::out | std::ofstream::trunc),
  logBuffer(logFile.rdbuf(), std::cout.rdbuf()),
  logLevel(logLevel_),
  indent(indent_),
  startTime(Time::getCurrentRealTime())
{
}

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
    fraction << std::fixed << double(time.getFractions()) / Time::FRACTIONS;
    std::cout << time.getSeconds() << fraction.str().substr(1,4) << " ";
    // next puts should go to logFile and std::out, done by this->put
    log = this;
  }
  if (category.length() > 0) {
    (*log) << category << ": ";
  } else {
    (*log) << sourceFileName.substr(9) << ':' << sourceFileLine << ": ";
  }
  return *log;
}

int Log::rank(-1);
std::string Log::fileName("cc4s.log");
unsigned int Log::logLevel(0);
LogStream *Log::logStream(nullptr);

void Log::setRank(int const rank_) {
  rank = rank_;
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

void Log::setLogLevel(const unsigned int logLevel_) {
  logLevel = logLevel_;
}

unsigned int Log::getLogLevel() {
  return logLevel;
}

LogStream &Log::getLogStream() {
  if (!logStream) {
    logStream = new LogStream(fileName, logLevel);
  }
  return *logStream;
}

