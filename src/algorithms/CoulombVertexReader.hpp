/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef COULOMB_VERTEX_READER_DEFINED
#define COULOMB_VERTEX_READER_DEFINED

#include <algorithms/Algorithm.hpp>
#include <ctf.hpp>
#include <cstdint>
#include <fstream>

namespace cc4s {
  class CoulombVertexReader: public Algorithm {
  public:
    ALGORITHM_REGISTRAR_DECLARATION(CoulombVertexReader);
    CoulombVertexReader(
      std::vector<Argument> const &argumentList
    );
    virtual ~CoulombVertexReader();
    virtual void run();
    virtual void dryRun();

  protected:
    int no, nv, nG, np;

    void readGammaGpqChunkBlocked(std::ifstream &file, CTF::Tensor<> &GammaGpq);
    void readEpsChunk(
      std::ifstream &file,  CTF::Tensor<> &epsi, CTF::Tensor<> &epsa
    );

    class Header {
    public:
      char magic[8];
      int32_t no, nv, nG, nSpins, kPoints, reserved_;
      static char const *MAGIC;
    };
    class Chunk {
    public:
      char magic[8];
      int64_t size;
      static char const *REALS_MAGIC;
      static char const *IMAGS_MAGIC;
      static char const *REALSIA_MAGIC;
      static char const *IMAGSIA_MAGIC;
      static char const *EPSILONS_MAGIC;
    };
  };
}

#endif

