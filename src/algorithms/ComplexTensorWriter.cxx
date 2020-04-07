/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <algorithms/ComplexTensorWriter.hpp>
#include <util/TensorIo.hpp>
#include <util/Log.hpp>
#include <fstream> 
#include <iomanip>
#include <ctf.hpp>
#include <util/Emitter.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(ComplexTensorWriter);

ComplexTensorWriter::ComplexTensorWriter(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

ComplexTensorWriter::~ComplexTensorWriter() {
}

void ComplexTensorWriter::run() {
  Tensor<complex> *A(getTensorArgument<complex>("Data"));
  std::string dataName(getArgumentData("Data")->getName());
  A->set_name(dataName.c_str());
  EMIT() << YAML::Key << "Data" << YAML::Value << dataName;
  std::string mode(getTextArgument("mode", "text"));
  if (mode == "binary") {
    // write binary
    std::string fileName(getTextArgument("file", dataName + ".bin"));
    TensorIo::writeBinary<complex>(fileName, *A);
    EMIT() << YAML::Key << "file" << YAML::Value << fileName;
  } else {
    // write text
    std::string fileName(getTextArgument("file", dataName + ".dat"));
    std::string rowIndexOrder(getTextArgument("rowIndexOrder", ""));
    std::string columnIndexOrder(getTextArgument("columnIndexOrder", ""));
    std::string delimiter(getTextArgument("delimiter", " "));
    TensorIo::writeText<complex>(
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

