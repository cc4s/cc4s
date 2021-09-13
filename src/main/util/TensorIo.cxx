#include <util/TensorIo.hpp>
#include <Reader.hpp>

// using namespace cc4s;

int cc4s::TensorIo::WRITE_REGISTERED =
  cc4s::Writer::registerWriteFunction("tensor", cc4s::TensorIo::write);

int cc4s::TensorIo::READ_REGISTERED =
  cc4s::Reader::registerReadFunction("tensor", cc4s::TensorIo::read);

