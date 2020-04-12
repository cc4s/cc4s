/*Copyright (c) 2018, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <algorithms/RealTensorReader.hpp>
#include <util/TensorIo.hpp>
#include <util/Emitter.hpp>
#include <util/Log.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>
#include <Data.hpp>

#include <fstream>

using namespace tcc;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(RealTensorReader);

RealTensorReader::RealTensorReader(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

void RealTensorReader::run() {
  std::string name(getArgumentData("Data")->getName());

  // make sure all processes start reading the file at the same time in case
  // it has been modified before
  Cc4s::world->barrier();
  setTensorArgument<Real<64>,DefaultTensorEngine>(
    "Data", read<Real<64>,DefaultTensorEngine>(name)
  );
}

template <typename F, typename TE>
PTR(ESC(tcc::Tensor<F,TE>)) RealTensorReader::read(const std::string &name) {
  PTR(ESC(tcc::Tensor<F,TE>)) A;
  std::string mode(getTextArgument("mode", "text"));
  if (mode == "binary") {
    std::string fileName(getTextArgument("file", name + ".bin"));
    EMIT() << YAML::Key << "file" << YAML::Value << fileName;
    A = TensorIo::readBinary<F,TE>(fileName);
  } else {
    std::string fileName(getTextArgument("file", name + ".dat").c_str());
    std::string delimiter(getTextArgument("delimiter", " "));
    int64_t bufferSize(getIntegerArgument("bufferSize", 128l*1024*1024));
    A = TensorIo::readText<F,TE>(fileName, delimiter, bufferSize);
    EMIT() << YAML::Key << "file" << YAML::Value << fileName;
  }
//  A->set_name(name.c_str());
  EMIT() << YAML::Key << "Data"  << YAML::Value << name;
  EMIT() << YAML::Key << "elements" << YAML::Value << A->getElementsCount();

  return A;
}

