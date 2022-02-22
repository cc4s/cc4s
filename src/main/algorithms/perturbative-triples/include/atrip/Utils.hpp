// [[file:~/cc4s/src/atrip/bbbfb30/atrip.org::*Prolog][Prolog:1]]
#pragma once
#include <sstream>
#include <string>
#include <map>
#include <chrono>

#include <ctf.hpp>
#include <atrip/Debug.hpp>

namespace atrip {
// Prolog:1 ends here

// [[file:~/cc4s/src/atrip/bbbfb30/atrip.org::*Pretty%20printing][Pretty printing:1]]
template <typename T>
  std::string pretty_print(T&& value) {
    std::stringstream stream;
#if ATRIP_DEBUG > 2
    dbg::pretty_print(stream, std::forward<T>(value));
#endif
    return stream.str();
  }
// Pretty printing:1 ends here

// [[file:~/cc4s/src/atrip/bbbfb30/atrip.org::*Chrono][Chrono:1]]
#define WITH_CHRONO(__chrono_name, ...)         \
  Atrip::chrono[__chrono_name].start();         \
  __VA_ARGS__                                   \
  Atrip::chrono[__chrono_name].stop();

struct Timer {
  using Clock = std::chrono::high_resolution_clock;
  using Event = std::chrono::time_point<Clock>;
  std::chrono::duration<double> duration;
  Event _start;
  inline void start() noexcept { _start = Clock::now(); }
  inline void stop() noexcept { duration += Clock::now() - _start; }
  inline void clear() noexcept { duration *= 0; }
  inline double count() const noexcept { return duration.count(); }
};
using Timings = std::map<std::string, Timer>;
// Chrono:1 ends here

// [[file:~/cc4s/src/atrip/bbbfb30/atrip.org::*Epilog][Epilog:1]]
}
// Epilog:1 ends here
