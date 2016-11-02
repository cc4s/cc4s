/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <algorithms/TensorWriter.hpp>
#include <util/TensorIo.hpp>
#include <util/Log.hpp>
#include <fstream> 
#include <iomanip>
#include <ctf.hpp>

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
  Tensor<> *A(getTensorArgument<>("Data"));
  std::string dataName(getArgumentData("Data")->getName());
  A->set_name(dataName.c_str());

  std::string mode(getTextArgument("mode", "text"));
  if (mode == "binary") {
    // write binary
    std::string fileName(getTextArgument("file", dataName + ".bin"));
    TensorIo::writeBinary<>(fileName, *A);
  } else {
    // write text
    std::string fileName(getTextArgument("file", dataName + ".dat"));
    std::string rowIndexOrder(getTextArgument("rowIndexOrder", ""));
    std::string columnIndexOrder(getTextArgument("columnIndexOrder", ""));
    std::string delimiter(getTextArgument("delimiter", " "));
    TensorIo::writeText<>(
      fileName, *A, rowIndexOrder, columnIndexOrder, delimiter
    );
  }
}

