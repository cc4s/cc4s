#ifndef TIME_DEFINED
#define TIME_DEFINED

#include <ctime>
#include <cstdint>
#include <ostream>
#include <iomanip>
#include <sstream>
#include <cmath>

namespace cc4s {
  /**
   * The Time class manages time intervals and their manipulation in terms
   * of seconds and nanoseconds. It is a wrapper of the timespec class of
   * ctime.
   */
  class Time: protected timespec {
  public:
    static constexpr int64_t FRACTIONS = 1000000000;
    static constexpr int FRACTION_DIGITS = 9;

    Time() {
      tv_sec = 0;
      tv_nsec = 0;
    }
    Time(int64_t seconds, int64_t nanoSeconds) {
      tv_sec = seconds;
      tv_nsec = nanoSeconds;
    }
    Time(Time const &t): timespec(t) {
    }

    int64_t getSeconds() const {
      return tv_sec;
    }
    int64_t getFractions() const {
      return tv_nsec;
    }
    double getFractionalSeconds() const {
      return tv_sec + static_cast<double>(tv_nsec) / FRACTIONS;
    }

    Time &operator += (Time const &t) {
      tv_nsec += tv_nsec;
      if (tv_nsec <= FRACTIONS) {
        tv_sec += t.tv_sec;
      } else {
        tv_sec += t.tv_sec + 1;
        tv_nsec -= FRACTIONS;
      }
      return *this;
    }

    Time &operator -= (Time const &t) {
      if (tv_nsec >= t.tv_nsec) {
        tv_nsec -= t.tv_nsec;
        tv_sec -= t.tv_sec;
      } else {
        tv_nsec += FRACTIONS - t.tv_nsec;
        tv_sec -= t.tv_sec + 1;
      }
      return *this;
    }

    static Time getCurrentRealTime() {
      Time time;
      clock_gettime(CLOCK_REALTIME, &time);
      return time;
    }
  };

  inline Time operator +(Time const &t1, Time const &t2) {
    Time result(t1);
    result += t2;
    return result;
  }

  inline Time operator -(Time const &t1, Time const &t2) {
    Time result(t1);
    result -= t2;
    return result;
  }

  inline Time operator *(const Time &difference, const double factor) {
    double d(difference.getFractionalSeconds() * factor);
    return Time(std::floor(d), std::round((d-std::floor(d)) * Time::FRACTIONS));
  }

  inline std::ostream &operator <<(std::ostream &stream, Time const &t) {
    std::stringstream time;
    time << t.getSeconds() << "."
      << std::setw(Time::FRACTION_DIGITS) << std::setfill('0')
      << t.getFractions();
    return stream << time.str();
  }
}

#endif

