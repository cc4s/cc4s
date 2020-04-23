#include <algorithms/SliceCoulombVertex.hpp>

#include <math/Complex.hpp>
#include <math/MathFunctions.hpp>
#include <util/Log.hpp>
#include <util/Exception.hpp>
#include <tcc/Tcc.hpp>
#include <Cc4s.hpp>

using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(SliceCoulombVertex);

Ptr<MapNode> SliceCoulombVertex::run(const Ptr<MapNode> &arguments) {
  // multiplex calls to template methods
  if (Cc4s::options->dryRun) {
    return run<DryTensorEngine>(arguments);
  } else {
    return run<DefaultTensorEngine>(arguments);
  }
}

template <typename TE>
Ptr<MapNode> SliceCoulombVertex::run(
  const Ptr<MapNode> &arguments
) {
  typedef Tensor<Complex<>,TE> T;
  auto coulombVertex(arguments->getMap("coulombVertex"));
  auto GammaGqr(coulombVertex->getValue<Ptr<T>>("data"));
  Assert(GammaGqr, "expecting coulombVertex to be complex");

  // read dimensions of hole and particle energies from meta data
  auto holeEnergies(arguments->getMap("holeEigenEnergies"));
  auto No(holeEnergies->getMap("dimensions")->getValue<size_t>(0));
  auto particleEnergies(arguments->getMap("particleEigenEnergies"));
  auto Nv(particleEnergies->getMap("dimensions")->getValue<size_t>(0));

  auto NG(GammaGqr->lens[0]);
  auto Np(GammaGqr->lens[1]);

  LOG(1,getName()) << "NG=" << NG << std::endl;
  LOG(1,getName()) << "No=" << No << std::endl;
  LOG(1,getName()) << "Nv=" << Nv << std::endl;
  LOG(1,getName()) << "Np=" << Np << std::endl;

  auto slices(New<MapNode>());
  // NOTE: only recipes returned, calculations done on demand
  {
    // NOTE: all result tensors assume shape on first assignment
    auto GammaGab(Tcc<TE>::template tensor<Complex<>>("GammaGab"));
    slices->setValue<Ptr<TensorRecipe<Complex<>,TE>>>(
      "Gab",
      COMPILE_RECIPE(GammaGab,
        (*GammaGab)["Gab"] <<= (*(*GammaGqr)({0,No,No},{NG,Np,Np}))["Gab"]
      )
    );
    // NOTE: result tensors, such as GammaGab, are kept inside recipes.
    // put in own scope { } to prevent accidental confusion with other recipes.
  }
  {
    auto GammaGai(Tcc<TE>::template tensor<Complex<>>("GammaGai"));
    slices->setValue<Ptr<TensorRecipe<Complex<>,TE>>>(
      "Gai",
      COMPILE_RECIPE(GammaGai,
        (*GammaGai)["Gai"] <<= (*(*GammaGqr)({0,No,0},{NG,Np,No}))["Gai"]
      )
    );
  }
  {
    auto GammaGia(Tcc<TE>::template tensor<Complex<>>("GammaGia"));
    slices->setValue<Ptr<TensorRecipe<Complex<>,TE>>>(
      "Gia",
      COMPILE_RECIPE(GammaGia,
        (*GammaGia)["Gia"] <<= (*(*GammaGqr)({0,0,No},{NG,No,Np}))["Gia"]
      )
    );
  }
  {
    auto GammaGij(Tcc<TE>::template tensor<Complex<>>("GammaGij"));
    slices->setValue<Ptr<TensorRecipe<Complex<>,TE>>>(
      "Gij",
      COMPILE_RECIPE(GammaGij,
        (*GammaGij)["Gij"] <<= (*(*GammaGqr)({0,0,0},{NG,No,No}))["Gij"]
      )
    );
  }

  // create result
  auto slicedCoulombVertex(New<MapNode>());
  slicedCoulombVertex->get("unit") = coulombVertex->get("unit");
  slicedCoulombVertex->get("spins") = coulombVertex->get("spins");
  slicedCoulombVertex->get("orbitals") = coulombVertex->get("orbitals");
  slicedCoulombVertex->get("slices") = slices;
  auto result(New<MapNode>());
  result->get("slicedCoulombVertex") = slicedCoulombVertex;
  return result;
}

