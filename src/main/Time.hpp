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

#ifndef TIME_DEFINED
#define TIME_DEFINED

#include <Real.hpp>

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
    Real<> getFractionalSeconds() const {
      return tv_sec + static_cast<Real<>>(tv_nsec) / FRACTIONS;
    }

    Time &operator += (Time const &t) {
      tv_nsec += t.tv_nsec;
      if (tv_nsec < FRACTIONS) {
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

  inline Time operator *(const Time &difference, const Real<> factor) {
    Real<> d(difference.getFractionalSeconds() * factor);
    return Time(floor(d), round((d-floor(d)) * Time::FRACTIONS));
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

