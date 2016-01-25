/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef CC4S_DEFINED
#define CC4S_DEFINED

#include <ctf.hpp>
#include <Options.hpp>
#include <Chi.hpp>
#include <CoulombIntegrals.hpp>
#include <Amplitudes.hpp>
#include <sys/time.h>

namespace cc4s {
  class Cc4s {
  public:
    Cc4s();
    ~Cc4s();
    void run();

    // static properties, accessible from everywhere
    static CTF::World *world;
    static Options *options;
    static Chi *chiReal, *chiImag;
    static CoulombIntegrals *V;
    static Amplitudes *T;

  protected:
/*
    void mp2Algorithm();
    void rpaAlgorithm();
*/
    void iterateMp2();
    void iterateCcsd();
    void iterateRpa();
    void iterateRccd();
    void iterateRccsd();

    void printBanner();
    void printStatistics();
    // TODO: move in proper class
    double diff(timespec const &start, timespec const &end);

    CTF::Flop_counter flopCounter;
  };
}

#endif

