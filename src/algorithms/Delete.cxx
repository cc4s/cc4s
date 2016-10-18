/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <algorithms/Delete.hpp>
#include <ctf.hpp>

using namespace CTF;
using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(Delete);

Delete::Delete(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

Delete::~Delete() {
}

void Delete::run() {
  Data *data(getArgumentData("Data"));
  if (data) {
    std::string dataName(data->getName());
    delete data;
    // remention it in case it will be written to in the future
    new Data(dataName);
  } else {
    LOG(0, "Delete") << "Data not allocated." << std::endl;
  }
}

void Delete::dryRun() {
  Data *data(getArgumentData("Data"));
  if (data) {
    std::string dataName(data->getName());
    delete data;
    // remention it in case it will be written to in the future
    new Data(dataName);
  } else {
    LOG(0, "Delete") << "Data not allocated." << std::endl;
  }
}

