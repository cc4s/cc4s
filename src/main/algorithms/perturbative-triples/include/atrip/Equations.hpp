// [[file:~/cc4s/src/atrip/bbbfb30/atrip.org::*Equations][Equations:1]]
#pragma once

#include<atrip/Slice.hpp>
#include<atrip/Blas.hpp>

namespace atrip {

  template <typename F=double>
  double getEnergyDistinct
    ( const F epsabc
    , std::vector<F> const& epsi
    , std::vector<F> const& Tijk_
    , std::vector<F> const& Zijk_
    ) {
    constexpr size_t blockSize=16;
    F energy(0.);
    const size_t No = epsi.size();
    for (size_t kk=0; kk<No; kk+=blockSize){
      const size_t kend( std::min(No, kk+blockSize) );
      for (size_t jj(kk); jj<No; jj+=blockSize){
        const size_t jend( std::min( No, jj+blockSize) );
        for (size_t ii(jj); ii<No; ii+=blockSize){
          const size_t iend( std::min( No, ii+blockSize) );
          for (size_t k(kk); k < kend; k++){
            const F ek(epsi[k]);
            const size_t jstart = jj > k ? jj : k;
            for (size_t j(jstart); j < jend; j++){
              F const ej(epsi[j]);
              F const facjk = j == k ? F(0.5) : F(1.0);
              size_t istart = ii > j ? ii : j;
              for (size_t i(istart); i < iend; i++){
                const F
                    ei(epsi[i])
                  , facij = i == j ? F(0.5) : F(1.0)
                  , denominator(epsabc - ei - ej - ek)
                  , U(Zijk_[i + No*j + No*No*k])
                  , V(Zijk_[i + No*k + No*No*j])
                  , W(Zijk_[j + No*i + No*No*k])
                  , X(Zijk_[j + No*k + No*No*i])
                  , Y(Zijk_[k + No*i + No*No*j])
                  , Z(Zijk_[k + No*j + No*No*i])
                  , A(maybeConjugate<F>(Tijk_[i + No*j + No*No*k]))
                  , B(maybeConjugate<F>(Tijk_[i + No*k + No*No*j]))
                  , C(maybeConjugate<F>(Tijk_[j + No*i + No*No*k]))
                  , D(maybeConjugate<F>(Tijk_[j + No*k + No*No*i]))
                  , E(maybeConjugate<F>(Tijk_[k + No*i + No*No*j]))
                  , F(maybeConjugate<F>(Tijk_[k + No*j + No*No*i]))
                  , value
                    = 3.0 * ( A * U
                              + B * V
                              + C * W
                              + D * X
                              + E * Y
                              + F * Z )
                   + ( ( U + X + Y )
                     - 2.0 * ( V + W + Z )
                     ) * ( A + D + E )
                   + ( ( V + W + Z )
                     - 2.0 * ( U + X + Y )
                     ) * ( B + C + F )
                  ;
                energy += 2.0 * value / denominator * facjk * facij;
              } // i
            } // j
          } // k
        } // ii
      } // jj
    } // kk
    return std::real(energy);
  }


  template <typename F=double>
  double getEnergySame
    ( const F epsabc
    , std::vector<F> const& epsi
    , std::vector<F> const& Tijk_
    , std::vector<F> const& Zijk_
    ) {
    constexpr size_t blockSize = 16;
    const size_t No = epsi.size();
    F energy = F(0.);
    for (size_t kk=0; kk<No; kk+=blockSize){
      const size_t kend( std::min( kk+blockSize, No) );
      for (size_t jj(kk); jj<No; jj+=blockSize){
        const size_t jend( std::min( jj+blockSize, No) );
        for (size_t ii(jj); ii<No; ii+=blockSize){
          const size_t iend( std::min( ii+blockSize, No) );
          for (size_t k(kk); k < kend; k++){
            const F ek(epsi[k]);
            const size_t jstart = jj > k ? jj : k;
            for(size_t j(jstart); j < jend; j++){
              const F facjk( j == k ? F(0.5) : F(1.0));
              const F ej(epsi[j]);
              const size_t istart = ii > j ? ii : j;
              for(size_t i(istart); i < iend; i++){
                const F
                  ei(epsi[i])
                , facij ( i==j ? F(0.5) : F(1.0))
                , denominator(epsabc - ei - ej - ek)
                , U(Zijk_[i + No*j + No*No*k])
                , V(Zijk_[j + No*k + No*No*i])
                , W(Zijk_[k + No*i + No*No*j])
                , A(maybeConjugate<F>(Tijk_[i + No*j + No*No*k]))
                , B(maybeConjugate<F>(Tijk_[j + No*k + No*No*i]))
                , C(maybeConjugate<F>(Tijk_[k + No*i + No*No*j]))
                , value
                  = F(3.0) * ( A * U
                             + B * V
                             + C * W
                             )
                  - ( A + B + C ) * ( U + V + W )
                ;
                energy += F(2.0) * value / denominator * facjk * facij;
              } // i
            } // j
          } // k
        } // ii
      } // jj
    } // kk
    return std::real(energy);
  }

  template <typename F=double>
  void singlesContribution
    ( size_t No
    , size_t Nv
    , const ABCTuple &abc
    , F const* Tph
    , F const* VABij
    , F const* VACij
    , F const* VBCij
    , F *Zijk
    ) {
    const size_t a(abc[0]), b(abc[1]), c(abc[2]);
    for (size_t k=0; k < No; k++)
    for (size_t i=0; i < No; i++)
    for (size_t j=0; j < No; j++) {
      const size_t ijk = i + j*No + k*No*No
                ,  jk = j + No * k
                ;
      Zijk[ijk] += Tph[ a + i * Nv ] * VBCij[ j + k * No ];
      Zijk[ijk] += Tph[ b + j * Nv ] * VACij[ i + k * No ];
      Zijk[ijk] += Tph[ c + k * Nv ] * VABij[ i + j * No ];
    }
  }

  template <typename F=double>
  void doublesContribution
    ( const ABCTuple &abc
    , size_t const No
    , size_t const Nv
    // -- VABCI
    , F const* VABph
    , F const* VACph
    , F const* VBCph
    , F const* VBAph
    , F const* VCAph
    , F const* VCBph
    // -- VHHHA
    , F const* VhhhA
    , F const* VhhhB
    , F const* VhhhC
    // -- TA
    , F const* TAphh
    , F const* TBphh
    , F const* TCphh
    // -- TABIJ
    , F const* TABhh
    , F const* TAChh
    , F const* TBChh
    // -- TIJK
    , F *Tijk
    ) {

    const size_t a = abc[0], b = abc[1], c = abc[2]
              , NoNo = No*No, NoNv = No*Nv
              ;

  #if defined(ATRIP_USE_DGEMM)
  #define _IJK_(i, j, k) i + j*No + k*NoNo
  #define REORDER(__II, __JJ, __KK)                                 \
    WITH_CHRONO("doubles:reorder",                                  \
    for (size_t k = 0; k < No; k++)                                 \
    for (size_t j = 0; j < No; j++)                                 \
    for (size_t i = 0; i < No; i++) {                               \
      Tijk[_IJK_(i, j, k)] += _t_buffer[_IJK_(__II, __JJ, __KK)];   \
    }                                                               \
    )
  #define DGEMM_PARTICLES(__A, __B)      \
    atrip::xgemm<F>( "T"                 \
                   , "N"                 \
                   , (int const*)&NoNo   \
                   , (int const*)&No     \
                   , (int const*)&Nv     \
                   , &one                \
                   , __A                 \
                   , (int const*)&Nv     \
                   , __B                 \
                   , (int const*)&Nv     \
                   , &zero               \
                   , _t_buffer.data()    \
                   , (int const*)&NoNo   \
                   );
  #define DGEMM_HOLES(__A, __B, __TRANSB)    \
    atrip::xgemm<F>( "N"                     \
                   , __TRANSB                \
                   , (int const*)&NoNo       \
                   , (int const*)&No         \
                   , (int const*)&No         \
                   , &m_one                  \
                   , __A                     \
                   , (int const*)&NoNo       \
                   , __B                     \
                   , (int const*)&No         \
                   , &zero                   \
                   , _t_buffer.data()        \
                   , (int const*)&NoNo       \
                   );
  #define MAYBE_CONJ(_conj, _buffer)                 \
    for (size_t __i = 0; __i < NoNoNo; ++__i)        \
      _conj[__i] = maybeConjugate<F>(_buffer[__i]);  \

    const size_t NoNoNo = No*NoNo;
    std::vector<F> _t_buffer;
    _t_buffer.reserve(NoNoNo);
    F one{1.0}, m_one{-1.0}, zero{0.0};

    WITH_CHRONO("double:reorder",
      for (size_t k = 0; k < NoNoNo; k++) {
         Tijk[k] = 0.0;
       })

    // TOMERGE: replace chronos
    WITH_CHRONO("doubles:holes",
      { // Holes part ========================================================

        std::vector<F> _vhhh(NoNoNo);

        // VhhhC[i + k*No + L*NoNo] * TABhh[L + j*No]; H1
        MAYBE_CONJ(_vhhh, VhhhC)
        WITH_CHRONO("doubles:holes:1",
          DGEMM_HOLES(_vhhh.data(), TABhh, "N")
          REORDER(i, k, j)
        )
        // VhhhC[j + k*No + L*NoNo] * TABhh[i + L*No]; H0
        WITH_CHRONO("doubles:holes:2",
          DGEMM_HOLES(_vhhh.data(), TABhh, "T")
          REORDER(j, k, i)
        )

        // VhhhB[i + j*No + L*NoNo] * TAChh[L + k*No]; H5
        MAYBE_CONJ(_vhhh, VhhhB)
        WITH_CHRONO("doubles:holes:3",
          DGEMM_HOLES(_vhhh.data(), TAChh, "N")
          REORDER(i, j, k)
        )
        // VhhhB[k + j*No + L*NoNo] * TAChh[i + L*No]; H3
        WITH_CHRONO("doubles:holes:4",
          DGEMM_HOLES(_vhhh.data(), TAChh, "T")
          REORDER(k, j, i)
        )

        // VhhhA[j + i*No + L*NoNo] * TBChh[L + k*No]; H1
        MAYBE_CONJ(_vhhh, VhhhA)
        WITH_CHRONO("doubles:holes:5",
          DGEMM_HOLES(_vhhh.data(), TBChh, "N")
          REORDER(j, i, k)
        )
        // VhhhA[k + i*No + L*NoNo] * TBChh[j + L*No]; H4
        WITH_CHRONO("doubles:holes:6",
          DGEMM_HOLES(_vhhh.data(), TBChh, "T")
          REORDER(k, i, j)
        )

      }
    )
  #undef MAYBE_CONJ

    WITH_CHRONO("doubles:particles",
      { // Particle part =====================================================
        // TAphh[E + i*Nv + j*NoNv] * VBCph[E + k*Nv]; P0
        WITH_CHRONO("doubles:particles:1",
          DGEMM_PARTICLES(TAphh, VBCph)
          REORDER(i, j, k)
        )
        // TAphh[E + i*Nv + k*NoNv] * VCBph[E + j*Nv]; P3
        WITH_CHRONO("doubles:particles:2",
          DGEMM_PARTICLES(TAphh, VCBph)
          REORDER(i, k, j)
        )
        // TCphh[E + k*Nv + i*NoNv] * VABph[E + j*Nv]; P5
        WITH_CHRONO("doubles:particles:3",
          DGEMM_PARTICLES(TCphh, VABph)
          REORDER(k, i, j)
        )
        // TCphh[E + k*Nv + j*NoNv] * VBAph[E + i*Nv]; P2
        WITH_CHRONO("doubles:particles:4",
          DGEMM_PARTICLES(TCphh, VBAph)
          REORDER(k, j, i)
        )
        // TBphh[E + j*Nv + i*NoNv] * VACph[E + k*Nv]; P1
        WITH_CHRONO("doubles:particles:5",
          DGEMM_PARTICLES(TBphh, VACph)
          REORDER(j, i, k)
        )
        // TBphh[E + j*Nv + k*NoNv] * VCAph[E + i*Nv]; P4
        WITH_CHRONO("doubles:particles:6",
          DGEMM_PARTICLES(TBphh, VCAph)
          REORDER(j, k, i)
        )
      }
    )

  #undef REORDER
  #undef DGEMM_HOLES
  #undef DGEMM_PARTICLES
  #undef _IJK_
  #else
    for (size_t k = 0; k < No; k++)
    for (size_t j = 0; j < No; j++)
    for (size_t i = 0; i < No; i++){
      const size_t ijk = i + j*No + k*NoNo
                ,  jk = j + k*No
                ;
      Tijk[ijk] = 0.0; // :important
      // HOLE DIAGRAMS: TABHH and VHHHA
      for (size_t L = 0; L < No; L++){
        // t[abLj] * V[Lcik]        H1
        // t[baLi] * V[Lcjk]        H0      TODO: conjugate T for complex
        Tijk[ijk] -= TABhh[L + j*No] * VhhhC[i + k*No + L*NoNo];
        Tijk[ijk] -= TABhh[i + L*No] * VhhhC[j + k*No + L*NoNo];

        // t[acLk] * V[Lbij]        H5
        // t[caLi] * V[Lbkj]        H3
        Tijk[ijk] -= TAChh[L + k*No] * VhhhB[i + j*No + L*NoNo];
        Tijk[ijk] -= TAChh[i + L*No] * VhhhB[k + j*No + L*NoNo];

        // t[bcLk] * V[Laji]        H2
        // t[cbLj] * V[Laki]        H4
        Tijk[ijk] -= TBChh[L + k*No] * VhhhA[j + i*No + L*NoNo];
        Tijk[ijk] -= TBChh[j + L*No] * VhhhA[k + i*No + L*NoNo];
      }
      // PARTILCE DIAGRAMS: TAPHH and VABPH
      for (size_t E = 0; E < Nv; E++) {
        // t[aEij] * V[bcEk]        P0
        // t[aEik] * V[cbEj]        P3 // TODO: CHECK THIS ONE, I DONT KNOW
        Tijk[ijk] += TAphh[E + i*Nv + j*NoNv] * VBCph[E + k*Nv];
        Tijk[ijk] += TAphh[E + i*Nv + k*NoNv] * VCBph[E + j*Nv];

        // t[cEki] * V[abEj]        P5
        // t[cEkj] * V[baEi]        P2
        Tijk[ijk] += TCphh[E + k*Nv + i*NoNv] * VABph[E + j*Nv];
        Tijk[ijk] += TCphh[E + k*Nv + j*NoNv] * VBAph[E + i*Nv];

        // t[bEji] * V[acEk]        P1
        // t[bEjk] * V[caEi]        P4
        Tijk[ijk] += TBphh[E + j*Nv + i*NoNv] * VACph[E + k*Nv];
        Tijk[ijk] += TBphh[E + j*Nv + k*NoNv] * VCAph[E + i*Nv];
      }

    }
  #endif
  }

}
// Equations:1 ends here
