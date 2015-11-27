#ifndef LOG_DEFINED
#define LOG_DEFINED

#include <string>
#include <iostream>
#include <ctf.hpp>

namespace cc4s {
  /**
   * \brief Class with static members offering control over logging.
   * Log entries are created with the macro LOG.
   */
  class Log {
  public:
    /**
     * \brief The log level to use for subsequent LOG messages.
     * A log message will only be
     * written if its log level is equal or below the current log level.
     */
    static int logLevel;
    /**
     * \brief Indentation string used for each log level.
     * By default a tab character will be used.
     */
    static std::string indent;

    static std::ostream &logStream(int const level);
  };

  template <typename F>
  void logMatrix(int level, CTF::Matrix<F> &matrix);
}

/**
 * \brief Provides an output stream for writing a log message of the
 * log level specified by the argument.
 * Note that this macro must be used as a statement and cannot be used as an
 * rvalue.
 */
#define LOG(level) \
  if (level > cc4s::Log::logLevel) { } else cc4s::Log::logStream(level)

#endif

