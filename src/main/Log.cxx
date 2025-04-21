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

#include <Complex.hpp>
#include <Log.hpp>

#include <sstream>
#include <string>

using namespace cc4s;

// Initialize static Log members
int Log::rank(-1);
size_t Log::flush_log_after_seconds = 0;
std::string Log::fileName("cc4s.log");
std::ofstream Log::logStream;
std::ofstream Log::nullStream;
std::stringstream Log::log_buffer;
Log::Duration Log::log_time_since_last_flush = Log::Duration(0);
Log::Event Log::log_last_now = Log::Clock::now();
Log::HeaderFunction
    /* header functions */
    Log::outHeaderFunction = [](const SourceLocation &) { return ""; },
    Log::errorHeaderFunction =
        [](const SourceLocation &location) {
          std::stringstream header;
          header << location << ": ERROR: ";
          return header.str();
        },
    Log::warningHeaderFunction =
        [](const SourceLocation &location) {
          std::stringstream header;
          header << location << ": WARNING: ";
          return header.str();
        },
    Log::logHeaderFunction = [](const SourceLocation &) { return ""; };

// Implementation of the member functions
Log::HeaderFunction Log::getLogHeaderFunction() noexcept {
  return logHeaderFunction;
}

void Log::setLogHeaderFunction(const Log::HeaderFunction &f) noexcept {
  logHeaderFunction = f;
}

Log::HeaderFunction Log::getWarningHeaderFunction() noexcept {
  return warningHeaderFunction;
}

void Log::setWarningHeaderFunction(const Log::HeaderFunction &f) noexcept {
  warningHeaderFunction = f;
}

Log::HeaderFunction Log::getErrorHeaderFunction() noexcept {
  return errorHeaderFunction;
}

void Log::setErrorHeaderFunction(const Log::HeaderFunction &f) noexcept {
  errorHeaderFunction = f;
}

Log::HeaderFunction Log::getOutHeaderFunction() noexcept {
  return outHeaderFunction;
}

void Log::setOutHeaderFunction(const Log::HeaderFunction &f) noexcept {
  outHeaderFunction = f;
}

void Log::flushLogBuffer() noexcept {
  // dump the buff to logStream
  // note for future: rdbuf also empties the buffer
  logStream << log_buffer.rdbuf();
}

std::ostream &Log::getLogStream() noexcept {
  if (rank == 0) {
    log_time_since_last_flush = Clock::now() - log_last_now;
    if (std::chrono::duration_cast<std::chrono::seconds>(
            log_time_since_last_flush) >=
        std::chrono::seconds(flush_log_after_seconds)) {
      // set the duration to 0
      log_time_since_last_flush = Duration(0);
      // set last flush to now
      log_last_now = Clock::now();
      flushLogBuffer();
    }
    return log_buffer;
  } else {
    return nullStream;
  }
}

std::ostream &Log::getOutStream() noexcept {
  if (rank == 0 && !Cc4s::dryRun)
    return std::cout;
  else
    return nullStream;
}

std::ostream &Log::getErrorStream() noexcept {
  if (rank == 0)
    return std::cerr;
  else
    return nullStream;
}

std::ostream &Log::getWarningStream() noexcept {
  if (rank == 0)
    return std::cout;
  else
    return nullStream;
}

// initialize logging
void Log::setRank(int const rank_) noexcept {
  rank = rank_;
  if (rank == 0) {
    logStream.open(fileName.c_str(), std::ofstream::out | std::ofstream::trunc);
  }
  // keep a nul stream should the output be masked, e.g. on non-root ranks
  nullStream.setstate(std::ios_base::badbit);
}

void Log::closeLogFile() noexcept {
  if (rank == 0) {
    log_buffer << "//;; Closing log file\n";
    flushLogBuffer();
    logStream.close();
  }
}

int Log::getRank() noexcept { return rank; }

void Log::setFileName(const std::string &fileName_) noexcept {
  fileName = fileName_;
}

std::string Log::getFileName() noexcept { return fileName; }
