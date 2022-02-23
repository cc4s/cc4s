// [[file:~/cc4s/src/atrip/bbbfb30/atrip.org::*Unions][Unions:1]]
#pragma once
#include <atrip/SliceUnion.hpp>

namespace atrip {

  template <typename F=double>
  void sliceIntoVector
    ( std::vector<F> &v
    , CTF::Tensor<F> &toSlice
    , std::vector<int64_t> const low
    , std::vector<int64_t> const up
    , CTF::Tensor<F> const& origin
    , std::vector<int64_t> const originLow
    , std::vector<int64_t> const originUp
    ) {
    // Thank you CTF for forcing me to do this
    struct { std::vector<int> up, low; }
        toSlice_ = { {up.begin(), up.end()}
                   , {low.begin(), low.end()} }
      , origin_ = { {originUp.begin(), originUp.end()}
                  , {originLow.begin(), originLow.end()} }
      ;

    WITH_OCD
    WITH_RANK << "slicing into " << pretty_print(toSlice_.up)
                          << "," << pretty_print(toSlice_.low)
              << " from " << pretty_print(origin_.up)
                   << "," << pretty_print(origin_.low)
              << "\n";

#ifndef ATRIP_DONT_SLICE
    toSlice.slice( toSlice_.low.data()
                 , toSlice_.up.data()
                 , 0.0
                 , origin
                 , origin_.low.data()
                 , origin_.up.data()
                 , 1.0);
    memcpy(v.data(), toSlice.data, sizeof(F) * v.size());
#endif

  }


  template <typename F=double>
  struct TAPHH : public SliceUnion<F> {
    TAPHH( CTF::Tensor<F> const& sourceTensor
         , size_t No
         , size_t Nv
         , size_t np
         , MPI_Comm child_world
         , MPI_Comm global_world
         ) : SliceUnion<F>( sourceTensor
                          , {Slice<F>::A, Slice<F>::B, Slice<F>::C}
                          , {Nv, No, No} // size of the slices
                          , {Nv}
                          , np
                          , child_world
                          , global_world
                          , Slice<F>::TA
                          , 6) {
           this->init(sourceTensor);
         }

    void sliceIntoBuffer(size_t it, CTF::Tensor<F> &to, CTF::Tensor<F> const& from) override
    {
      const int Nv = this->sliceLength[0]
              , No = this->sliceLength[1]
              , a = this->rankMap.find({static_cast<size_t>(Atrip::rank), it});
              ;


      sliceIntoVector<F>( this->sources[it]
                        , to,   {0, 0, 0},    {Nv, No, No}
                        , from, {a, 0, 0, 0}, {a+1, Nv, No, No}
                        );

    }

  };


  template <typename F=double>
  struct HHHA : public SliceUnion<F> {
    HHHA( CTF::Tensor<F> const& sourceTensor
        , size_t No
        , size_t Nv
        , size_t np
        , MPI_Comm child_world
        , MPI_Comm global_world
        ) : SliceUnion<F>( sourceTensor
                         , {Slice<F>::A, Slice<F>::B, Slice<F>::C}
                         , {No, No, No} // size of the slices
                         , {Nv}         // size of the parametrization
                         , np
                         , child_world
                         , global_world
                         , Slice<F>::VIJKA
                         , 6) {
           this->init(sourceTensor);
         }

    void sliceIntoBuffer(size_t it, CTF::Tensor<F> &to, CTF::Tensor<F> const& from) override
    {

      const int No = this->sliceLength[0]
              , a = this->rankMap.find({static_cast<size_t>(Atrip::rank), it})
              ;

      sliceIntoVector<F>( this->sources[it]
                        , to,   {0, 0, 0},    {No, No, No}
                        , from, {0, 0, 0, a}, {No, No, No, a+1}
                        );

    }
  };

  template <typename F=double>
  struct ABPH : public SliceUnion<F> {
    ABPH( CTF::Tensor<F> const& sourceTensor
        , size_t No
        , size_t Nv
        , size_t np
        , MPI_Comm child_world
        , MPI_Comm global_world
        ) : SliceUnion<F>( sourceTensor
                         , { Slice<F>::AB, Slice<F>::BC, Slice<F>::AC
                           , Slice<F>::BA, Slice<F>::CB, Slice<F>::CA
                           }
                         , {Nv, No} // size of the slices
                         , {Nv, Nv} // size of the parametrization
                         , np
                         , child_world
                         , global_world
                         , Slice<F>::VABCI
                         , 2*6) {
           this->init(sourceTensor);
         }

    void sliceIntoBuffer(size_t it, CTF::Tensor<F> &to, CTF::Tensor<F> const& from) override {

      const int Nv = this->sliceLength[0]
              , No = this->sliceLength[1]
              , el = this->rankMap.find({static_cast<size_t>(Atrip::rank), it})
              , a = el % Nv
              , b = el / Nv
              ;


      sliceIntoVector<F>( this->sources[it]
                        , to,   {0, 0},       {Nv, No}
                        , from, {a, b, 0, 0}, {a+1, b+1, Nv, No}
                        );

    }

  };

  template <typename F=double>
  struct ABHH : public SliceUnion<F> {
    ABHH( CTF::Tensor<F> const& sourceTensor
        , size_t No
        , size_t Nv
        , size_t np
        , MPI_Comm child_world
        , MPI_Comm global_world
        ) : SliceUnion<F>( sourceTensor
                         , {Slice<F>::AB, Slice<F>::BC, Slice<F>::AC}
                         , {No, No} // size of the slices
                         , {Nv, Nv} // size of the parametrization
                         , np
                         , child_world
                         , global_world
                         , Slice<F>::VABIJ
                         , 6) {
           this->init(sourceTensor);
         }

    void sliceIntoBuffer(size_t it, CTF::Tensor<F> &to, CTF::Tensor<F> const& from) override {

      const int Nv = from.lens[0]
              , No = this->sliceLength[1]
              , el = this->rankMap.find({static_cast<size_t>(Atrip::rank), it})
              , a = el % Nv
              , b = el / Nv
              ;

      sliceIntoVector<F>( this->sources[it]
                        , to,   {0, 0},       {No, No}
                        , from, {a, b, 0, 0}, {a+1, b+1, No, No}
                        );


    }

  };


  template <typename F=double>
  struct TABHH : public SliceUnion<F> {
    TABHH( CTF::Tensor<F> const& sourceTensor
         , size_t No
         , size_t Nv
         , size_t np
         , MPI_Comm child_world
         , MPI_Comm global_world
         ) : SliceUnion<F>( sourceTensor
                          , {Slice<F>::AB, Slice<F>::BC, Slice<F>::AC}
                          , {No, No} // size of the slices
                          , {Nv, Nv} // size of the parametrization
                          , np
                          , child_world
                          , global_world
                          , Slice<F>::TABIJ
                          , 6) {
           this->init(sourceTensor);
         }

    void sliceIntoBuffer(size_t it, CTF::Tensor<F> &to, CTF::Tensor<F> const& from) override {
      // TODO: maybe generalize this with ABHH

      const int Nv = from.lens[0]
              , No = this->sliceLength[1]
              , el = this->rankMap.find({static_cast<size_t>(Atrip::rank), it})
              , a = el % Nv
              , b = el / Nv
              ;

      sliceIntoVector<F>( this->sources[it]
                        , to,   {0, 0},       {No, No}
                        , from, {a, b, 0, 0}, {a+1, b+1, No, No}
                        );


    }

  };

}
// Unions:1 ends here
