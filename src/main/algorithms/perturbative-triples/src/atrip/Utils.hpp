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

#include <atrip/Debug.hpp>

namespace atrip {
// Prolog:1 ends here

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
