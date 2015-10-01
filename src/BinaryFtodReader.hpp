/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef BINARY_FTOD_READER_DEFINED
#define BINARY_FTOD_READER_DEFINED

#include "FtodReader.hpp"
#include <cstdint>

class BinaryFtodReader: public FtodReader {
  public:
    virtual void read();
    virtual void write();

  protected:
    class Header {
       char magic[8];
       int32_t no, nv, nG, nSpins, kPoints, reserved_;
    };
    class Chunk {
      char magic[8];
      int64_t size;
    };
};


#endif
