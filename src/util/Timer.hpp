#ifndef TIMER_DEFINED
#define TIMER_DEFINED

#include <ctime>

namespace cc4s {
  /**
   * Timer class prividing timing functionality.
   * The object can be created in a scope whose lifetime
   * is then measured.
   */
  class Timer {
  public:
    Timer(double *seconds);
    ~Timer();
  protected:
    double *seconds;
    timespec start;
  };
}

#endif

