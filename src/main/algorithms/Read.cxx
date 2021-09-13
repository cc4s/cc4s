#include <algorithms/Read.hpp>

#include <Cc4s.hpp>
#include <Reader.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(Read)

/**
 * \brief Interface to Reader.
 */
Ptr<MapNode> Read::run(const Ptr<MapNode> &arguments) {
  auto fileName(arguments->getValue<std::string>("fileName"));

  try {
    auto data(Reader(fileName).read());
    // create result
    auto result(New<MapNode>(data->sourceLocation));
    result->get("data") = data;
    return result;
  } catch (const Ptr<Exception> &cause) {
    throw New<Exception>(
      "Failed to read file '" + fileName + "'",
      arguments->get("fileName")->sourceLocation,
      cause
    );
  }
}

