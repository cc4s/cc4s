/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#ifndef FTOD_READER_DEFINED
#define FTOD_READER_DEFINED

namespace cc4s {
  /**
   * \deprecated This will no longer be used in the Data-Algorithm
   * design.
   */
  class FtodReader {
  public:
    virtual void read() = 0;
  };
}

#endif

