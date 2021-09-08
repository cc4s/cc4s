#ifndef BLACS_WORLD_DEFINED
#define BLACS_WORLD_DEFINED

namespace cc4s {
  class BlacsWorld {
  public:
    BlacsWorld(int rank, int processes, int processRows = -1);
    ~BlacsWorld();
    void barrier();

    int rank;
    int context;
    int lens[2], firstElement[2];
  };
}

#endif

