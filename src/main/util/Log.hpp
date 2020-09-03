/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef LOG_DEFINED
#define LOG_DEFINED

#include <util/SourceLocation.hpp>

#include <string>
#include <iostream>
#include <streambuf>
#include <fstream>
#include <functional>

namespace cc4s {
  /**
   * \brief Class with static members offering control over logging.
   * Log entries are created with the macro LOG.
   */
  class Log {
  public:
    typedef std::function<std::string(const SourceLocation &)> HeaderFunction;

    static void setRank(const int rank);
    static int getRank();
    static void setFileName(const std::string &fileName);
    static std::string getFileName();
    static std::ofstream &getStream() {
      return stream;
    }
    static void setOutHeaderFunction(const HeaderFunction &f) {
      outHeaderFunction = f;
    }
    static HeaderFunction getOutHeaderFunction() {
      return outHeaderFunction;
    }
    static void setLogHeaderFunction(const HeaderFunction &f) {
      logHeaderFunction = f;
    }
    static HeaderFunction getLogHeaderFunction() {
      return logHeaderFunction;
    }

  protected:
    static int rank;
    static std::string fileName;
    static std::ofstream stream;
    static HeaderFunction outHeaderFunction, logHeaderFunction;
  };
}

/**
 * \brief Provides an output stream for writing a log message.
 */
#define OUT() \
  (std::cout << Log::getOutHeaderFunction()(SOURCE_LOCATION))
#define OUT_LOCATION(LOCATION) \
  (std::cout << Log::getOutHeaderFunction()(LOCATION))
#define LOG() \
  (Log::getStream() << Log::getLogHeaderFunction()(SOURCE_LOCATION))
#define LOG_LOCATION(LOCATION) \
  (Log::getStream() << Log::getLogHeaderFunction()(LOCATION))

#endif

