/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
    static std::ostream &getOutStream() {
      if (rank == 0) return std::cout;
      else return logStream;
    }
    static void setFileName(const std::string &fileName);
    static std::string getFileName();
    static std::ofstream &getLogStream() {
      return logStream;
    }
    static void setOutHeaderFunction(const HeaderFunction &f) {
      outHeaderFunction = f;
    }
    static HeaderFunction getOutHeaderFunction() {
      return outHeaderFunction;
    }
    static void setErrorHeaderFunction(const HeaderFunction &f) {
      errorHeaderFunction = f;
    }
    static HeaderFunction getErrorHeaderFunction() {
      return errorHeaderFunction;
    }
    static void setWarningHeaderFunction(const HeaderFunction &f) {
      warningHeaderFunction = f;
    }
    static HeaderFunction getWarningHeaderFunction() {
      return warningHeaderFunction;
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
    static std::ofstream logStream;
    static HeaderFunction
      outHeaderFunction, errorHeaderFunction, warningHeaderFunction,
      logHeaderFunction;
  };
}

/**
 * \brief Provides output streams for writing various types of messages
 * to the console and for writing to the log file.
 */
#define OUT() \
  (Log::getOutStream() << Log::getOutHeaderFunction()(SOURCE_LOCATION))
#define OUT_LOCATION(LOCATION) \
  (Log::getOutStream() << Log::getOutHeaderFunction()(LOCATION))
#define ERROR() \
  (Log::getOutStream() << Log::getErrorHeaderFunction()(SOURCE_LOCATION))
#define ERROR_LOCATION(LOCATION) \
  (Log::getOutStream() << Log::getErrorHeaderFunction()(LOCATION))
#define WARNING() \
  (Log::getOutStream() << Log::getWarningHeaderFunction()(SOURCE_LOCATION))
#define WARNING_LOCATION(LOCATION) \
  (Log::getOutStream() << Log::getWarningHeaderFunction()(LOCATION))
#define LOG() \
  (Log::getLogStream() << Log::getLogHeaderFunction()(SOURCE_LOCATION))
#define LOG_LOCATION(LOCATION) \
  (Log::getLogStream() << Log::getLogHeaderFunction()(LOCATION))

#endif

