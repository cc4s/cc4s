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

// [[file:~/cuda/atrip/atrip.org::*Prolog][Prolog:1]]
#pragma once
#include <sstream>
#include <string>

#define IGNORABLE(param)                                                       \
  do { (void)(param); } while (0)

#if defined(HAVE_CUDA)
#  include <nvToolsExt.h>
#endif /* defined(HAVE_CUDA) */

#if defined(__NVCC__)
#  pragma nv_diagnostic_push
#  if defined __NVCC_DIAG_PRAGMA_SUPPORT__
// http://www.ssl.berkeley.edu/~jimm/grizzly_docs/SSL/opt/intel/cc/9.0/lib/locale/en_US/mcpcom.msg
#    pragma nv_diag_suppress partial_override
#  else
#    pragma diag_suppress partial_override
#  endif
#  include <ctf.hpp>
#  pragma nv_diagnostic_pop
#else
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wvla"
#  pragma GCC diagnostic ignored "-Wnonnull"
#  pragma GCC diagnostic ignored "-Wall"
#  pragma GCC diagnostic ignored "-Wint-in-bool-context"
#  pragma GCC diagnostic ignored "-Wunused-parameter"
#  pragma GCC diagnostic ignored "-Wdeprecated-copy"
#  include <ctf.hpp>
#  pragma GCC diagnostic pop
#endif

#include <atrip/Debug.hpp>

namespace atrip {
// Prolog:1 ends here

// [[file:~/cuda/atrip/atrip.org::*Pretty%20printing][Pretty printing:1]]
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
template <typename T>
std::string pretty_print(T &&value) {
#if ATRIP_DEBUG > 2
  std::stringstream stream;
  dbg::pretty_print(stream, std::forward<T>(value));
  return stream.str();
#else
  return "";
#endif
}
#pragma GCC diagnostic pop
// Pretty printing:1 ends here

//  A nice handy macro to do formatting
#define _FORMAT(_fmt, ...)                                                     \
  ([&](void) -> std::string {                                                  \
    int _sz = std::snprintf(nullptr, 0, _fmt, __VA_ARGS__);                    \
    std::vector<char> _out(_sz + 1);                                           \
    std::snprintf(&_out[0], _out.size(), _fmt, __VA_ARGS__);                   \
    return std::string(_out.data());                                           \
  })()

// [[file:~/cuda/atrip/atrip.org::*Epilog][Epilog:1]]
} // namespace atrip
// Epilog:1 ends here
