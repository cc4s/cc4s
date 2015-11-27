/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef BINARY_FTOD_READER_DEFINED
#define BINARY_FTOD_READER_DEFINED

#include <FtodReader.hpp>
#include <Chi.hpp>
#include <cstdint>
#include <fstream>

namespace cc4s {
  class BinaryFtodReader: public FtodReader {
    public:
      virtual void read();
      virtual void write();

    protected:
      int no, nv, nG;
      int64_t np;
      void readChiChunk(std::ifstream &file, Chi *chi);
      void readEpsChunk(std::ifstream &file);

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
          static char const *EPSILONS_MAGIC;
      };
  };
}


#endif
