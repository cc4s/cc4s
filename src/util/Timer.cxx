#include <util/Timer.hpp>

#include <ctime>

using namespace cc4s;

/**
 * \brief Constructor of a Timer object starting the measurement.
 * \parameter seconds specifies where the measured time shall be
 * written to upon destruction of the object marking the end of the
 * measurement.
 */
Timer::Timer(double *seconds_): seconds(seconds_) {
  clock_gettime(CLOCK_REALTIME, &start);
}

/**
 * Destroys the timer object concluding the time measurement
 * and writing the elapsed time in the double specified upon
 * construction of the timer object.
 */
Timer::~Timer() {
  timespec end;
  clock_gettime(CLOCK_REALTIME, &end);
  *seconds = end.tv_sec-start.tv_sec + (end.tv_nsec-start.tv_nsec) / 1e9;
}
