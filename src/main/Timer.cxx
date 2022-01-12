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

#include <Timer.hpp>

using namespace cc4s;

/**
 * \brief Constructor of a Timer object starting the measurement.
 * \parameter seconds specifies where the measured time shall be
 * written to upon destruction of the object marking the end of the
 * measurement.
 */
Timer::Timer(Time *time_): time(time_), start(Time::getCurrentRealTime()) {
}

/**
 * Destroys the timer object concluding the time measurement
 * and adding the elapsed time to the time specified upon
 * construction of the timer object.
 */
Timer::~Timer() {
  Time end(Time::getCurrentRealTime());
  *time += end - start;
}
