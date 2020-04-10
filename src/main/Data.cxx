/*Copyright (c) 2015, Andreas Grueneis and Felix Hummel, all rights reserved.*/

#include <Data.hpp>
#include <sstream>

using namespace cc4s;

Data::Data(std::string const &name_): name(name_), stage(MENTIONED) {
  std::stringstream sStream;
  sStream << name_ << " of yet unknown type";
  typeName = sStream.str();
  dataMap[name_] = this;
}


std::map<std::string, Data *> Data::dataMap;

int TypedData::nextId;

