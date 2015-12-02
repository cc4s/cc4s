/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef PARTICLE_HOLE_COULOMB_VERTEX_READER_DEFINED
#define PARTICLE_HOLE_COULOMB_VERTEX_READER_DEFINED

#include <Algorithm.hpp>
#include <ctf.hpp>
#include <cstdint>
#include <fstream>

namespace cc4s {
  class ParticleHoleCoulombVertexReader: public Algorithm {
    public:
      ParticleHoleCoulombVertexReader(
        std::vector<Argument const *> const &argumentList
      );
      virtual ~ParticleHoleCoulombVertexReader();
      virtual void run();

    protected:
      int no, nv, nG;
      int64_t np;

      void readChiAiChunkBlocked(std::ifstream &file, CTF::Tensor<> *chiAi);
      void readEpsChunk(std::ifstream &file,  CTF::Tensor<> *ieps, CTF::Tensor<> *aeps);

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
