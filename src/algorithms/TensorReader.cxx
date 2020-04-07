/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <algorithms/TensorReader.hpp>
#include <util/TensorIo.hpp>
#include <util/Log.hpp>
#include <Cc4s.hpp>
#include <fstream>
#include <ctf.hpp>
#include <util/Emitter.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(TensorReader);

TensorReader::TensorReader(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

TensorReader::~TensorReader() {
}

void TensorReader::run() {
  std::string name(getArgumentData("Data")->getName());

  // make sure all processes start reading the file at the same time in case
  // it has been modified before
  MPI_Barrier(Cc4s::world->comm);

  int64_t precision(getIntegerArgument("precision", 64));
  switch (precision) {
  case 64:
    allocatedTensorArgument<Float64>(
      "Data", read<Float64>(name)
    );
    break;
  case 128:
#ifndef INTEL_COMPILER
    allocatedTensorArgument<Float128>(
      "Data", read<Float128>(name)
    );
#else
    throw new EXCEPTION("Quadruple precision not supported for Intel");
#endif
    break;
  }
}

template <typename F>
Tensor<F> *TensorReader::read(const std::string &name) {
  Tensor<F> *A;
  std::string mode(getTextArgument("mode", "text"));
  if (mode == "binary") {
    std::string fileName(getTextArgument("file", name + ".bin"));
    EMIT() << YAML::Key << "file" << YAML::Value << fileName;
    A = TensorIo::readBinary<F>(fileName);
  } else {
    std::string fileName(getTextArgument("file", name + ".dat").c_str());
    std::string delimiter(getTextArgument("delimiter", " "));
    int64_t bufferSize(getIntegerArgument("bufferSize", 128l*1024*1024));
    A = TensorIo::readText<F>(fileName, delimiter, bufferSize);
    EMIT() << YAML::Key << "file" << YAML::Value << fileName;
  }
  A->set_name(name.c_str());
  EMIT() << YAML::Key << "Data"  << YAML::Value << name;

  int64_t indexCount(1);
  for (int dim(0); dim < A->order; ++dim) {
    indexCount *= A->lens[dim];
  }
  EMIT() << YAML::Key << "elements" << YAML::Value << indexCount;

  return A;
}

