// [[file:~/cc4s/src/atrip/bbbfb30/atrip.org::*Header][Header:1]]
#pragma once
#include <sstream>
#include <string>
#include <map>
#include <chrono>

#include <ctf.hpp>

#include <atrip/Utils.hpp>

#define ADD_ATTRIBUTE(_type, _name, _default)   \
  _type _name = _default;                       \
  Input& with_ ## _name(_type i) {              \
    _name = i;                                  \
    return *this;                               \
  }

namespace atrip {

  struct Atrip {

    static int rank;
    static int np;
    static Timings chrono;
    static void init();

    template <typename F=double>
    struct Input {
      CTF::Tensor<F> *ei = nullptr
                        , *ea = nullptr
                        , *Tph = nullptr
                        , *Tpphh = nullptr
                        , *Vpphh = nullptr
                        , *Vhhhp = nullptr
                        , *Vppph = nullptr
                        ;
      Input& with_epsilon_i(CTF::Tensor<F> * t) { ei = t; return *this; }
      Input& with_epsilon_a(CTF::Tensor<F> * t) { ea = t; return *this; }
      Input& with_Tai(CTF::Tensor<F> * t) { Tph = t; return *this; }
      Input& with_Tabij(CTF::Tensor<F> * t) { Tpphh = t; return *this; }
      Input& with_Vabij(CTF::Tensor<F> * t) { Vpphh = t; return *this; }
      Input& with_Vijka(CTF::Tensor<F> * t) { Vhhhp = t; return *this; }
      Input& with_Vabci(CTF::Tensor<F> * t) { Vppph = t; return *this; }

      enum TuplesDistribution {
        NAIVE,
        GROUP_AND_SORT,
      };

      ADD_ATTRIBUTE(bool, rankRoundRobin, false)
      ADD_ATTRIBUTE(bool, chrono, false)
      ADD_ATTRIBUTE(bool, barrier, false)
      ADD_ATTRIBUTE(int, maxIterations, 0)
      ADD_ATTRIBUTE(int, iterationMod, -1)
      ADD_ATTRIBUTE(int, percentageMod, -1)
      ADD_ATTRIBUTE(TuplesDistribution, tuplesDistribution, NAIVE)

    };

    struct Output {
      double energy;
    };
    template <typename F=double>
    static Output run(Input<F> const& in);
  };

}
// Header:1 ends here
