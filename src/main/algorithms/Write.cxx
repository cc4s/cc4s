/*Copyright (c) 2020, Andreas Grueneis and Felix Hummel, all rights reserved.*/
#include <algorithms/Write.hpp>

#include <Cc4s.hpp>
#include <Writer.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(Write)

/**
 * \brief Interface to Writer.
 */
Ptr<MapNode> Write::run(const Ptr<MapNode> &arguments) {
  auto data(arguments->get("data"));
  ASSERT_LOCATION(data, "expecting key 'data'", arguments->sourceLocation);
  auto fileName(arguments->getValue<std::string>("fileName"));
  auto useBinary(arguments->getValue<bool>("binary", false));

  auto persistentData(Writer(fileName, useBinary).write(data));

  // create result
  auto result(New<MapNode>(data->sourceLocation));
  result->get("persitentData") = persistentData;
  return result;
}

