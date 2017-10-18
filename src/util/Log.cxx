/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include "Log.hpp"

#include <math/Complex.hpp>
#include <sstream>
#include <string>

using namespace cc4s;

LogStream::LogStream(
  std::string const &logFileName,
  int const logLevel_,
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
  int const rank,
  std::string const &sourceFileName,
  int const level,
  std::string const &category
) {
  Time time(Time::getCurrentRealTime());
  time -= startTime;
  logFile << time << " ";
  std::stringstream timeStr;
  timeStr << time.getSeconds() << "."
          << std::setw(3) << std::setfill('0')
          << std::to_string(time.getFractions()).substr(0,3);
  *this << timeStr.str() << " ";
  std::ostream *log(&logFile);
  if (logLevel >= level) {
    for (int i(0); i < level; ++i) {
      std::cout << indent.c_str();
    }
    // next puts should go to logFile and std::out, done by this->put
    log = this;
  }
  if (category == "root") logFile << "root: ";
  else (*log) << (category.length() > 0 ? category : sourceFileName) << ": ";
  return *log;
}

int Log::rank(-1);
LogStream *Log::logStream(nullptr);

void Log::setRank(int const rank_) {
  rank = rank_;
}

int Log::getRank() {
  return rank;
}

void Log::setLogStream(LogStream *logStream_) {
  logStream = logStream_;
}

LogStream &Log::getLogStream() {
  return *logStream;
}

