#ifndef LOG_DEFINED
#define LOG_DEFINED

#include <ctf.hpp>

namespace cc4s {
  /**
   * \brief Class with static members offering control over logging.
   * Log entries are created with the macro LOG.
   */
  class Log {
  public:
    /**
     * The log level to use for subsequent LOG messages.
     * A log message will only be
     * written if its log level is equal or below the current log level.
     */
    static int logLevel;
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
#define LOG(level) if (level > cc4s::Log::logLevel) { } else std::cout

#endif

