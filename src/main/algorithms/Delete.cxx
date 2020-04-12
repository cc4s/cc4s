/*Copyright (c) 2016, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <algorithms/Delete.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(Delete);

Delete::Delete(
  std::vector<Argument> const &argumentList
): Algorithm(argumentList) {
}

void Delete::run() {
  auto data(getArgumentData("Data"));
  if (data) {
    std::string dataName(data->getName());
    // remention empty data in case it will be written to in the future
    // this will overwrite and discard any previous data
    NEW(Data,dataName);
  } else {
    LOG(0, "Delete") << "Data not allocated." << std::endl;
  }
}

void Delete::dryRun() {
  auto data(getArgumentData("Data"));
  if (data) {
    std::string dataName(data->getName());
    // remention empty data in case it will be written to in the future
    // this will overwrite and discard any previous data
    NEW(Data,dataName);
  } else {
    LOG(0, "Delete") << "Data not allocated." << std::endl;
  }
}

