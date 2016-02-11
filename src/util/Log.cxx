#include "Log.hpp"

#include <util/Complex.hpp>

using namespace cc4s;

LogStream::LogStream(
  std::string const &logFileName,
  int const logLevel_,
  std::string const &indent_
):
  std::ostream(&logBuffer),
  logFile(logFileName.c_str()), logLevel(logLevel_), indent(indent_),
  logBuffer(logFile.rdbuf(), std::cout.rdbuf())
{
}

std::ostream &LogStream::prepare(int rank, int const level) {
  // TODO: add time here
  logFile << rank << " " << level << ": ";
  if (logLevel >= level) {
    for (int i(0); i < level; ++i) {
      std::cout << indent.c_str();
    }
    // next puts should go to logFile and std::out, done by this->put
    return *this;
  } else {
    // next puts should only go to logFile
    return logFile;
  }
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
  if (logStream) delete logStream;
  logStream = logStream_;
}

LogStream &Log::getLogStream() {
  return *logStream;
}

template <typename F>
void cc4s::logMatrix(int level, CTF::Matrix<F> &m) {
  F *values(new F[m.lens[0]*m.lens[1]]);
  m.read_all(values);
  for (int i(0); i < m.lens[0]; ++i) {
    for (int j(0); j < m.lens[1]; ++j) {
      LOG(level) << " " << values[i+j*m.lens[0]];
    }
    LOG(level) << std::endl;
  }
}

// instantiate
template
void cc4s::logMatrix(int level, CTF::Matrix<double> &m);
template
void cc4s::logMatrix(int level, CTF::Matrix<complex> &m);

