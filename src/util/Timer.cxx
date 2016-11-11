#include <util/Timer.hpp>

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
 * and writing the elapsed time in the double specified upon
 * construction of the timer object.
 */
Timer::~Timer() {
  *time = Time::getCurrentRealTime();
  *time -= start;
}
