#ifndef CHRONO_HPP_
#define CHRONO_HPP_

#include <map>
#include <chrono>
#include <string>

// [[file:~/cuda/atrip/atrip.org::*Chrono][Chrono:1]]
#if defined(HAVE_CUDA)
#  define WITH_CHRONO(__chrono_name, ...)                                      \
    nvtxRangePushA(__chrono_name);                                             \
    Atrip::chrono[__chrono_name].start();                                      \
    __VA_ARGS__                                                                \
    Atrip::chrono[__chrono_name].stop();                                       \
    nvtxRangePop();
#else
#  define WITH_CHRONO(__chrono_name, ...)                                      \
    Atrip::chrono[__chrono_name].start();                                      \
    __VA_ARGS__                                                                \
    Atrip::chrono[__chrono_name].stop();
#endif

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

#endif
