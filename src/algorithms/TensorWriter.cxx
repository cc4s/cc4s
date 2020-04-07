/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <algorithms/TensorWriter.hpp>
#include <util/TensorIo.hpp>
#include <util/Log.hpp>
#include <fstream> 
#include <iomanip>
#include <ctf.hpp>
#include <util/Emitter.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(TensorWriter);

TensorWriter::TensorWriter(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

TensorWriter::~TensorWriter() {
}

void TensorWriter::run() {
  std::string dataName(getArgumentData("Data")->getName());

  int64_t precision(getIntegerArgument("precision", 64));
  switch (precision) {
  case 64:
    write<Float64>(dataName);
    break;
  case 128:
#ifndef INTEL_COMPILER
    write<Float128>(dataName);
#else
    throw new EXCEPTION("Quadruple precision not supported for Intel");
#endif
    break;
  }
}

template <typename F>
void TensorWriter::write(const std::string &name) {
  Tensor<F> *A(getTensorArgument<F>("Data"));
  A->set_name(name.c_str());
  EMIT() << YAML::Key << "Data" << YAML::Value << name;
  std::string mode(getTextArgument("mode", "text"));
  if (mode == "binary") {
    // write binary
    std::string fileName(getTextArgument("file", name + ".bin"));
    TensorIo::writeBinary<F>(fileName, *A);
    EMIT() << YAML::Key << "file" << YAML::Value << fileName;
  } else {
    // write text
    std::string fileName(getTextArgument("file", name + ".dat"));
    std::string rowIndexOrder(getTextArgument("rowIndexOrder", ""));
    std::string columnIndexOrder(getTextArgument("columnIndexOrder", ""));
    std::string delimiter(getTextArgument("delimiter", " "));
    TensorIo::writeText<F>(
      fileName, *A, rowIndexOrder, columnIndexOrder, delimiter
    );
    EMIT() << YAML::Key << "file" << YAML::Value << fileName;
  }

  int64_t indexCount(1);
  for (int dim(0); dim < A->order; ++dim) {
    indexCount *= A->lens[dim];
  }
  EMIT() << YAML::Key << "elements" << YAML::Value << indexCount;

}

