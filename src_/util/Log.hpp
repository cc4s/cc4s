/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LOG_DEFINED
#define LOG_DEFINED

#include <util/Time.hpp>

#include <string>
#include <iostream>
#include <streambuf>
#include <fstream>

namespace cc4s {
  class LogBuffer: public std::streambuf {
  public:
    LogBuffer(
      std::streambuf *log_, std::streambuf *out_
    ): log(log_), out(out_) { }
  protected:
    virtual int overflow(int c) {
      if (c == EOF) {
        return !EOF;
      } else {
        int const logPut(log->sputc(c));
        int const outPut(out->sputc(c));
        return logPut == EOF || outPut == EOF ? EOF : c;
      }
    }

    virtual int sync() {
      int const logSync(log->pubsync());
      int const outSync(out->pubsync());
      return logSync == 0 && outSync == 0 ? 0 : -1;
    }   
    std::streambuf *log, *out;
  };

  class LogStream: public std::ostream {
  public:
    LogStream(
      std::string const &logFileName,
      int const logLevel = 0,
      std::string const &indent = "\t"
    );

    std::ostream &prepare(
      int const rank,
      std::string const &sourceFileName,
      int const level,
      std::string const &category = ""
    );
  protected:
    std::ofstream logFile;
    LogBuffer logBuffer;

    /**
     * \brief The log level to use for subsequent LOG messages.
     * A log message will only be
     * written if its log level is equal or below the current log level.
     */
    int logLevel;
    /**
     * \brief Indentation string used for each log level.
     * By default a tab character will be used.
     */
    std::string indent;

    Time startTime;
  };

  /**
   * \brief Class with static members offering control over logging.
   * Log entries are created with the macro LOG.
   */
  class Log {
  public:
    static void setRank(const int rank);
    static int getRank();
    static void setFileName(const std::string &fileName);
    static std::string getFileName();
    static void setLogLevel(const int logLevel);
    static int getLogLevel();

    static LogStream &getLogStream();

  protected:
    static int rank;
    static std::string fileName;
    static int logLevel;
    static LogStream *logStream;
  };
}

// TODO: return output stream for all processes, including those
// who shouldn't print, so that formating functions can be used.
/**
 * \brief Provides an output stream for writing a log message of the
 * log level specified by the argument.
 * Note that this macro must be used as a statement and cannot be used as an
 * rvalue.
 */
#define OUT() \
  if (cc4s::Log::getRank() != 0) { \
  } else std::cout
#define NEW_FILE(NAME) \
  if (cc4s::Log::getRank() != 0) { \
  } else std::ofstream(NAME, std::ofstream::out)
#define FILE(NAME) \
  if (cc4s::Log::getRank() != 0) { \
  } else std::ofstream(NAME, std::ofstream::app)
#define LOG(...) \
  if (cc4s::Log::getRank() != 0) { \
  } else cc4s::Log::getLogStream().prepare(0, __FILE__, __VA_ARGS__)
#define LOG_RANK(...) \
  cc4s::Log::getLogStream().prepare(cc4s::Log::getRank(), __FILE__, __VA_ARGS__)

#endif

