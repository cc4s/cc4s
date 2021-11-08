#include <algorithms/DimensionProperty.hpp>

#include <Cc4s.hpp>
#include <Parser.hpp>
#include <util/Scanner.hpp>
#include <tcc/Tcc.hpp>
#include <math/Real.hpp>
#include <math/Complex.hpp>
#include <math/MathFunctions.hpp>
#include <vector>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(DimensionProperty)

/**
 * \brief Extracts a property matrix for a property of a
 * selected dimension of the given tensor.
 * The first dimension of the matrix corresponds
 * to the selected dimension of the given tensor,
 * the second dimension corresponds to each possible property value.
 * The matrix contains ones for all combinations where
 * the property of the index in the first dimension matches the
 * property value of the second dimension.
 * The matrix is zero otherwise.
 */
Ptr<MapNode> DimensionProperty::run(const Ptr<MapNode> &arguments) {
  auto op(arguments->getMap("operator"));
  auto scalarType(op->getValue<std::string>("scalarType"));
  // multiplex calls to template methods
  if (Cc4s::options->dryRun) {
    using TE = DefaultDryTensorEngine;
    if (scalarType == TypeTraits<Real<>>::getName()) {  
      return run<Real<>,TE>(arguments);
    } else if (scalarType == TypeTraits<Complex<>>::getName()) {
      return run<Complex<>,TE>(arguments);
    }
  } else {
    using TE = DefaultTensorEngine;
    if (scalarType == TypeTraits<Real<>>::getName()) {  
      return run<Real<>,TE>(arguments);
    } else if (scalarType == TypeTraits<Complex<>>::getName()) {
      return run<Complex<>,TE>(arguments);
    }
  }
  ASSERT_LOCATION(
    false, "scalar type '" + scalarType + "' not supported",
    op->get("scalarType")->sourceLocation
  );
}

template <typename F, typename TE>
Ptr<MapNode> DimensionProperty::run(
  const Ptr<MapNode> &arguments
) {
  typedef Tensor<F,TE> T;
  auto op(arguments->getMap("operator"));
  auto data(op->getValue<Ptr<T>>("data"));
  ASSERT_LOCATION(
    data, "expecting operator to be a tensor", op->sourceLocation
  );

  auto dimensionIndex(arguments->getValue<Natural<>>("dimension"));
  auto propertyName(arguments->getValue<std::string>("property"));

  ASSERT_LOCATION(
    dimensionIndex < data->getLens().size(),
    "Tensor "
      + data->getName()
      + " has no dimension #"
      + std::to_string(dimensionIndex),
    op->sourceLocation
  );

  auto propertyIt(
    data->dimensions[dimensionIndex]->properties.find(propertyName)
  );
  ASSERT_LOCATION(
    propertyIt != data->dimensions[dimensionIndex]->properties.end(),
    "Tensor "
      + data->getName()
      + " dimension "
      + std::to_string(dimensionIndex)
      + " does not have a property "
      + propertyName,
    op->sourceLocation
  );
  auto property(propertyIt->second);

  // create property matrix
  std::vector<Natural<>> lens(
    {data->lens[dimensionIndex], property->indicesOfProperty.size()}
  );
  std::vector<std::string> types(
    {data->dimensions[dimensionIndex]->name, propertyName}
  );
  auto P( Tcc<TE>::template tensor<F>(lens, propertyName) );

  // enter one at the respective positions
  std::vector<Natural<>> indices(0);
  std::vector<F> values(0);
  for (Natural<> p(0); p < property->indicesOfProperty.size(); ++p) {
    for (Natural<> index: property->indicesOfProperty[p]) {
      if (Cc4s::world->getRank() == 0) {
        indices.push_back(index + lens[0] * p);
        values.push_back(F(1));
      }
    }
  }
  P->write(indices.size(), indices.data(), values.data());

  // build meta data of Delta tensor
  auto DimensionPropertyNode( New<MapNode>(op->sourceLocation) );
  DimensionPropertyNode->setValue<>("version", std::string("v1.0"));
  DimensionPropertyNode->setValue<>("scalarType", TypeTraits<F>::getName());
  auto dimensionsNode( New<MapNode>(SOURCE_LOCATION) );
  for (Natural<> d(0); d < lens.size(); ++d) {
    auto dimensionNode( New<MapNode>(SOURCE_LOCATION) );
    dimensionNode->setValue<>("length", lens[d]);
    dimensionNode->setValue<>("type", types[d]);
    // NOTE: no type written
    dimensionsNode->get(d) = dimensionNode;
  }
  DimensionPropertyNode->get("dimensions") = dimensionsNode;
  DimensionPropertyNode->setValue("data", P);

  auto result(New<MapNode>(SOURCE_LOCATION));
  result->get("dimensionProperty") = DimensionPropertyNode;
  return result;
}

