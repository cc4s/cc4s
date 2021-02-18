/*Copyright (c) 2021, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <math/TensorUnion.hpp>
#include <Writer.hpp>

 using namespace cc4s;

int TensorUnionIo::WRITE_REGISTERED =
  Writer::registerWriteFunction("tensorUnion", TensorUnionIo::write);

//int TensorIo::READ_REGISTERED =
//  Reader::registerReadFunction("tensor", cc4s::TensorIo::read);

