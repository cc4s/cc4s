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

#ifndef TIMER_DEFINED
#define TIMER_DEFINED

#include <util/Time.hpp>

namespace cc4s {
  /**
   * Timer class providing timing functionality.
   * The object can be created in a scope whose lifetime
   * is then measured.
   */
  class Timer {
  public:
    Timer(Time *time);
    ~Timer();
  protected:
    Time *time;
    Time start;
  };
}

#endif

