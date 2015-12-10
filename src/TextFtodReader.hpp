/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef TEXT_FTOD_READER_DEFINED
#define TEXT_FTOD_READER_DEFINED

#include <FtodReader.hpp>

namespace cc4s {
  class TextFtodReader: public FtodReader {
  public:
    virtual void read();
  };
}

#endif

