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
  delete getArgumentData("Data");
}

void Delete::dryRun() {
  delete getArgumentData("Data");
}

