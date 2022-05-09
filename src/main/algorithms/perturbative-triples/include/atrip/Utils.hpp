// Copyright 2022 Alejandro Gallo
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// [[file:../../atrip.org::*Prolog][Prolog:1]]
#pragma once
#include <sstream>
#include <string>
#include <map>
#include <chrono>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvla"
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wdeprecated-copy"
#include <ctf.hpp>
#pragma GCC diagnostic pop

#include <atrip/Debug.hpp>

namespace atrip {
// Prolog:1 ends here

// [[file:../../atrip.org::*Pretty printing][Pretty printing:1]]
template <typename T>
  std::string pretty_print(T&& value) {
    std::stringstream stream;
#if ATRIP_DEBUG > 2
    dbg::pretty_print(stream, std::forward<T>(value));
#endif
    return stream.str();
  }
// Pretty printing:1 ends here

// [[file:../../atrip.org::*Chrono][Chrono:1]]
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

// [[file:../../atrip.org::*Epilog][Epilog:1]]
}
// Epilog:1 ends here
