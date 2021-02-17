/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <util/TensorIo.hpp>
#include <Writer.hpp>

using namespace cc4s;

int TensorIo::WRITE_REGISTERED =
  Writer::registerWriteFunction("tensor", TensorIo::write);

