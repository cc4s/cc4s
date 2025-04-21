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

#include <SourceLocation.hpp>
#include <Cc4s.hpp>

#include <string>
#include <iostream>
#include <streambuf>
#include <fstream>
#include <functional>
#include <chrono>

namespace cc4s {
  /**
   * \brief Class with static members offering control over logging.
   * Log entries are created with the macro LOG.
   */
  class Log {
  public:
    // type alias 'HeaderFunction' for a function that takes a
    // SourceLocation and returns a string.
    typedef std::function<std::string(const SourceLocation &)> HeaderFunction;

    using Clock = std::chrono::high_resolution_clock;
    using Event = std::chrono::time_point<Clock>;
    using Duration = std::chrono::duration<double>;

    static void setRank(const int rank) noexcept;
    static void closeLogFile() noexcept;
    static int getRank() noexcept;
    static std::ostream &getOutStream() noexcept;
    static std::ostream &getErrorStream() noexcept;
    static std::ostream &getWarningStream() noexcept;
    static void setFileName(const std::string &fileName) noexcept;
    static std::string getFileName() noexcept;
    static std::ostream &getLogStream() noexcept;
    static void setOutHeaderFunction(const HeaderFunction &f) noexcept;
    static HeaderFunction getOutHeaderFunction() noexcept;
    static void setErrorHeaderFunction(const HeaderFunction &f) noexcept;
    static HeaderFunction getErrorHeaderFunction() noexcept;
    static void setWarningHeaderFunction(const HeaderFunction &f) noexcept;
    static HeaderFunction getWarningHeaderFunction() noexcept;
    static void setLogHeaderFunction(const HeaderFunction &f) noexcept;
    static HeaderFunction getLogHeaderFunction() noexcept;
    static void flushLogBuffer() noexcept;
    static size_t flush_log_after_seconds;

  protected:
    static Duration log_time_since_last_flush;
    static Event log_last_now;
    static std::stringstream log_buffer;
    static int rank;
    static std::string fileName;
    static std::ofstream logStream;
    // nullStream discards output from non-rank 0 processes, ensuring
    // only rank 0 logs.
    static std::ofstream nullStream;
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
  (Log::getErrorStream() << Log::getErrorHeaderFunction()(SOURCE_LOCATION))
#define ERROR_LOCATION(LOCATION) \
  (Log::getErrorStream() << Log::getErrorHeaderFunction()(LOCATION))
#define WARNING() \
  (Log::getErrorStream() << Log::getWarningHeaderFunction()(SOURCE_LOCATION))
#define WARNING_LOCATION(LOCATION) \
  (Log::getErrorStream() << Log::getWarningHeaderFunction()(LOCATION))
#define LOG() \
  (Log::getLogStream() << Log::getLogHeaderFunction()(SOURCE_LOCATION))
#define LOG_LOCATION(LOCATION) \
  (Log::getLogStream() << Log::getLogHeaderFunction()(LOCATION))

#endif
