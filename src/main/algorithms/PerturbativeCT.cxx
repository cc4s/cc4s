#include <memory>

#include <algorithms/PerturbativeCT.hpp>
#include <TensorSet.hpp>


using namespace cc4s;

ALGORITHM_REGISTRAR_DEFINITION(PerturbativeCT)

Ptr<MapNode> PerturbativeCT::run(const Ptr<MapNode> &arguments) {
  // multiplex calls to template methods
  Ptr<MapNode> result;
  if (Cc4s::dryRun) {
    using TE = DefaultDryTensorEngine;
    (
      result = runReal<TE>(arguments)
    ) || (
      result = runComplex<TE>(arguments)
    );
  } else {
    using TE = DefaultTensorEngine;
    (
      result = runReal<TE>(arguments)
    ) || (
      result = runComplex<TE>(arguments)
    );
  }
  ASSERT_LOCATION(
    result, "unsupported tensor type as 'coulombIntegrals'",
    arguments->sourceLocation
  );
  return result;
}

template <typename TE>
Ptr<MapNode> PerturbativeCT::runReal(
  const Ptr<MapNode> &arguments
) {

  bool noRing(
    arguments->template getValue<bool>("noRing", false)
  );


  using TSr = TensorSet<Real<>,TE>;
  using TSc = TensorSet<Complex<>,TE>;
  auto amplitudes( arguments->getPtr<TSr>("amplitudes") );
  if (!amplitudes) return nullptr;

  auto Tph ( amplitudes->get("ph") );
  auto Tpphh( amplitudes->get("pphh") );

  auto coulombIntegrals( arguments->getPtr<TSr>("coulombIntegrals") );

  auto Vpppp( coulombIntegrals->get("pppp") );
  auto Vhhhh( coulombIntegrals->get("hhhh") );
  auto Vhhhp( coulombIntegrals->get("hhhp") );
  auto Vhhpp( coulombIntegrals->get("hhpp") );
  auto Vhphp( coulombIntegrals->get("hphp") );
  auto Vhhph( coulombIntegrals->get("hhph") );
  auto Vhpph( coulombIntegrals->get("hpph") );
  auto Vhphh( coulombIntegrals->get("hphh") );
  auto Vpphh( coulombIntegrals->get("pphh") );


  // if Vppph is a real tensor we store a backup of it in Bppph
  // and restore it at the end of the routine
  // If it is a recipe, we simply construct Vppph in this routine
  auto Bppph( Tcc<TE>::template tensor<Real<>>("Bppph") );
  auto Vppph( Tcc<TE>::template tensor<Real<>>("Vppph") );
  auto dressedGph( Tcc<TE>::template tensor<Complex<>>("dressedGph") );


  if (arguments->isGiven("coulombVertex")) {
    auto coulombVertex( arguments->getPtr<TSc>("coulombVertex") );
    auto GammaGpp( coulombVertex->get("pp") );
    auto GammaGph( coulombVertex->get("ph") );

    COMPILE(
//      (*conjGpp)["Gqr"] <<=
//        cc4s::map<Complex<>>(conj<Complex<>>, (*GammaGpp)["Grq"]),
      (*dressedGph)["Gck"] <<=
        (*GammaGpp)["Gcf"] * map<Complex<>>(toComplex<Real<>>, (*Tph)["fk"])
    )->execute();

    auto realGammaGpp( Tcc<TE>::template tensor<Real<>>("realGammaGpp") );
    auto imagGammaGpp( Tcc<TE>::template tensor<Real<>>("imagGammaGpp") );
    auto realGammaGph( Tcc<TE>::template tensor<Real<>>("realGammaGph") );
    auto imagGammaGph( Tcc<TE>::template tensor<Real<>>("imagGammaGph") );
    COMPILE(
      (*realGammaGpp)["Gab"] <<= map<Real<>>(real<Complex<>>, (*GammaGpp)["Gab"] ),
      (*imagGammaGpp)["Gab"] <<= map<Real<>>(imag<Complex<>>, (*GammaGpp)["Gab"] ),
      (*realGammaGph)["Gai"] <<= map<Real<>>(real<Complex<>>, (*GammaGph)["Gai"] ),
      (*imagGammaGph)["Gai"] <<= map<Real<>>(imag<Complex<>>, (*GammaGph)["Gai"] ),
      (*Vppph)["bcek"] <<= (*realGammaGpp)["Geb"] * (*realGammaGph)["Gck"],
      (*Vppph)["bcek"]  += (*imagGammaGpp)["Geb"] * (*imagGammaGph)["Gck"]
    )->execute();

  } else {
    // now coulombIntegrals is a tensor (not a recipe)
    auto whatever(coulombIntegrals->get("ppph"));
    auto recipe(dynamicPtrCast<TensorRecipe<Real<>,TE>>(whatever));
    if (recipe) {
      std::string msg("No CoulombVertex given but Integrals appear as Recipe. Abort!");
      ASSERT_LOCATION(false, msg, SOURCE_LOCATION);
    }
    COMPILE(
      (*Vppph)["bcek"] <<= (*whatever)["bcek"],
      (*Bppph)["bcek"] <<= (*Vppph)["bcek"]
    )->execute();
  }

  auto Xijka( Tcc<TE>::template tensor<Real<>>("Xijka") );
  auto Iia( Tcc<TE>::template tensor<Real<>>("Iia") );
  auto Ippph( Tcc<TE>::template tensor<Real<>>("Ippph") );
  auto dpImcjk( Tcc<TE>::template tensor<Real<>>("dpImcjk") );
  auto dpIjkmc( Tcc<TE>::template tensor<Real<>>("dpIjkmc") );

  std::shared_ptr<TSr> J = std::make_shared<TSr>();
  J->get("ppph") = Ippph;
  J->get("hhhp") = dpIjkmc;



  // Implementing the equations from Piecuch, Computer Physics Communications 149, 71 (2002)
  // Equations in Table 1: we will write the corresponding line number
  // Note that we need only half of the Equations for the intermediates Jppph and Jhhhp

  // We contract the ppXiajb && ppXiabj right away. This reduces the overall
  // memory high water mark
  {
   auto dpXiajb( Tcc<TE>::template tensor<Real<>>("dpXiajb") );
   auto dpXiabj( Tcc<TE>::template tensor<Real<>>("dpXiabj") );
   COMPILE(
     // Line 7: Xijka
     (*Xijka)["ijka"] <<= (*Vhhhp)["ijka"],
     (*Xijka)["ijka"] += (1.0) * (*Tph)["ek"] * (*Vpphh)["eaij"],

     // Line 11: ppXiajb
     (*dpXiajb)["iajb"] <<= (*Vhphp)["iajb"],
     (*dpXiajb)["iajb"] += (-1.0) * (*Vhhhp)["imjb"] * (*Tph)["am"],
     //(*dpXiajb)["iajb"] += ( 0.5) * (*Tph)["ej"] * (*Ippph)["beai"],

     // Line 15: ppXiabj
     //(*dpXiabj)["iabj"] <<= (*Vhpph)["iabj"],
     (*dpXiabj)["iabj"] <<= (*Vpphh)["baij"],
     (*dpXiabj)["iabj"] += (-1.0) * (*Vhhph)["imbj"] * (*Tph)["am"],
     //(*dpXiabj)["iabj"] += ( 0.5) * (*Ippph)["ebai"] * (*Tph)["ej"]
     // Line 19: 1st to 4th term
     (*dpImcjk)["iajk"] <<= (*Vhphh)["iajk"],
     (*dpImcjk)["iajk"] += (-1.0) * (*Vhhhh)["imjk"] * (*Tph)["am"],
     (*dpImcjk)["iajk"] += ( 1.0) * (*dpXiajb)["iaje"] * (*Tph)["ek"],
     (*dpImcjk)["iajk"] += ( 1.0) * (*Tph)["ej"] * (*dpXiabj)["iaek"],
     // Line 19: 5th - 9th term
     (noRing) ? (
       Tcc<TE>::sequence()
     ) : (
       (*dpImcjk)["iajk"] += ( 2.0) * (*Xijka)["imje"] * (*Tpphh)["aekm"],
       Tcc<TE>::sequence()
     ),
     (*dpImcjk)["iajk"] += (-1.0) * (*Xijka)["imje"] * (*Tpphh)["eakm"],
     (*dpImcjk)["iajk"] += (-1.0) * (*Xijka)["mije"] * (*Tpphh)["eamk"],
     (*dpImcjk)["iajk"] += (-1.0) * (*Tpphh)["aemj"] * (*Xijka)["mike"]

   )->execute();

  }
  {
  auto pXiajb( Tcc<TE>::template tensor<Real<>>("pXiajb") );
  auto pXiabj( Tcc<TE>::template tensor<Real<>>("pXiabj") );
  COMPILE(
    // Line 6: pXaibc (we store it in a ppph object for memory reasons)
    (*Ippph)["bcai"] <<= (*Vppph)["bcai"],
    // Line 6: pXaibc
    (*Ippph)["bcai"] += (-0.5) * (*Tph)["am"] * (*Vpphh)["bcmi"],
    // Line 1: Iia
    (*Iia)["ia"] <<= ( 2.0) * (*Vpphh)["aeim"] * (*Tph)["em"],
    (*Iia)["ia"]  += (-1.0) * (*Vpphh)["eaim"] * (*Tph)["em"],
    // (*Iia)["ia"] += Fia["ia"], // TODO: use non-canonical orbitals

    // Line 10: pXiajb
    (*pXiajb)["iajb"] <<= (*Vhphp)["iajb"],
    (*pXiajb)["iajb"] += (-0.5) * (*Vhhhp)["imjb"] * (*Tph)["am"],
    (*pXiajb)["iajb"] += ( 1.0) * (*Tph)["ej"] * (*Ippph)["beai"],

    // Line 14: pXiabj
    //(*pXiabj)["iabj"] <<= (*Vhpph)["iabj"],
    (*pXiabj)["iabj"] <<= (*Vpphh)["baij"],
    (*pXiabj)["iabj"] += (-0.5) * (*Vhhph)["imbj"] * (*Tph)["am"],
    (*pXiabj)["iabj"] += ( 1.0) * (*Ippph)["ebai"] * (*Tph)["ej"],

    // Line 5: Xaibc (store it in the same object to safe memory)
    (*Ippph)["bcai"] += (-0.5) * (*Tph)["am"] * (*Vpphh)["bcmi"],

     // Line 17: 3rd term - 5th term
    (*Vppph)["abci"] += (-1.0) * (*pXiajb)["maic"] * (*Tph)["bm"],
    (*Vppph)["abci"] += (-1.0) * (*Tph)["am"] * (*pXiabj)["mbci"],
    (*Vppph)["abci"] += (-1.0) * (*Iia)["mc"] * (*Tpphh)["abmi"]

  )->execute();
  }
  // Line 17: we use Vppph as scratch for ppIabci
  if (arguments->isGiven("coulombVertex")) {
    auto coulombVertex( arguments->getPtr<TSc>("coulombVertex") );
    auto GammaGpp( coulombVertex->get("pp") );
    auto GammaGph( coulombVertex->get("ph") );
    auto realGammaGpp( Tcc<TE>::template tensor<Real<>>("realGammaGpp") );
    auto imagGammaGpp( Tcc<TE>::template tensor<Real<>>("imagGammaGpp") );
    auto realDressedGph( Tcc<TE>::template tensor<Real<>>("realDressedGph") );
    auto imagDressedGph( Tcc<TE>::template tensor<Real<>>("imagDressedGph") );
    COMPILE(
      (*realGammaGpp)["Gab"] <<= map<Real<>>(real<Complex<>>, (*GammaGpp)["Gab"] ),
      (*imagGammaGpp)["Gab"] <<= map<Real<>>(imag<Complex<>>, (*GammaGpp)["Gab"] ),
      (*realDressedGph)["Gai"] <<= map<Real<>>(real<Complex<>>, (*dressedGph)["Gai"] ),
      (*imagDressedGph)["Gai"] <<= map<Real<>>(imag<Complex<>>, (*dressedGph)["Gai"] ),
      //(*dpIbcek)["bcek"] += (*realConjGpp)["Gbe"] * (*realDressedGph)["Gck"],
      //(*dpIbcek)["bcek"] += (-1.0) * (*imagConjGpp)["Gbe"] * (*imagDressedGph)["Gck"]
      (*Vppph)["bcek"] += (*realGammaGpp)["Geb"] * (*realDressedGph)["Gck"],
      (*Vppph)["bcek"] += (*imagGammaGpp)["Geb"] * (*imagDressedGph)["Gck"]
    )->execute();

  } else {
    COMPILE(
      //(*dpIbcek)["bcek"] += ( 1.0) * (*Vpppp)["bcef"] * (*Tph)["fk"]
      (*Vppph)["bcek"] += ( 1.0) * (*Vpppp)["bcef"] * (*Tph)["fk"]
    )->execute();
  }

  cc4s::enable_lowmem(true);
  COMPILE(
//    // Line 17: 6th term - 10th term
    (noRing) ? (
      Tcc<TE>::sequence()
    ) : (
      (*Vppph)["abci"] += ( 2.0) * (*Ippph)["ceam"] * (*Tpphh)["beim"],
      Tcc<TE>::sequence()
    ),
    (*Vppph)["abci"] += (-1.0) * (*Ippph)["ceam"] * (*Tpphh)["ebim"],
    (*Vppph)["abci"] += (-1.0) * (*Ippph)["ecam"] * (*Tpphh)["ebmi"],
    (*Vppph)["abci"] += (-1.0) * (*Tpphh)["eaim"] * (*Ippph)["ecbm"]
  )->execute();
  cc4s::enable_lowmem(false);
  COMPILE(
    (*Vppph)["abci"] += ( 1.0) * (*Tpphh)["bamn"] * (*Xijka)["mnic"],
    // Line 19: 9th term
    (*dpImcjk)["iajk"] += ( 1.0) * (*Tpphh)["efkj"] * (*Ippph)["efai"],

    //Line 19: Ippph contributuion from Line 19: 3rd + 4th term
    (*dpImcjk)["iajk"] += ( 0.5) * (*Tph)["fj"] * (*Ippph)["efai"] * (*Tph)["ek"],
    (*dpImcjk)["iajk"] += ( 0.5) * (*Tph)["ej"] * (*Ippph)["feai"] * (*Tph)["fk"],
    // For technical reasons we want the result in the form Vhhhp
    (*dpIjkmc)["jkmc"] <<= (*dpImcjk)["mcjk"],
    (*Ippph)["bcek"] <<= (*Vppph)["bcek"]
  )->execute();

  // If Vppph was given as a tensor we have to restore the original result
  if (!arguments->isGiven("coulombVertex")) {
    COMPILE( (*Vppph)["bcek"] <<= (*Bppph)["bcek"] )->execute();
  }

  if (noRing) OUT() << "Evaluating ct-intermediate without ring contribution" << std::endl;
  //TODO: give tensor sets meta data.
  //auto metaData(New<MapNode>(SOURCE_LOCATION));
  //metaData->setValue("noRing", noRing);
  //J->metaData = metaData;
  auto result( New<MapNode>(SOURCE_LOCATION));
  result->setPtr("cTIntermediate", J);
  return result;

}

// The complex code is not yet memory optimized.
template <typename TE>
Ptr<MapNode> PerturbativeCT::runComplex(
  const Ptr<MapNode> &arguments
) {


  bool noRing(
    arguments->template getValue<bool>("noRing", false)
  );

  auto amplitudes( arguments->getPtr<TensorSet<Complex<>,TE>>("amplitudes") );
  if (!amplitudes) return nullptr;

  auto Tph ( amplitudes->get("ph") );
  auto Tpphh( amplitudes->get("pphh") );

  auto coulombIntegrals( arguments->getPtr<TensorSet<Complex<>,TE>>("coulombIntegrals") );

  auto Vpppp( coulombIntegrals->get("pppp") );
  auto Vhhhh( coulombIntegrals->get("hhhh") );
  auto Vhhhp( coulombIntegrals->get("hhhp") );
  auto Vhhpp( coulombIntegrals->get("hhpp") );
  auto Vhphp( coulombIntegrals->get("hphp") );
  auto Vhhph( coulombIntegrals->get("hhph") );
  auto Vhpph( coulombIntegrals->get("hpph") );
  auto Vhphh( coulombIntegrals->get("hphh") );


  // if Vppph is a real tensor we store a backup of it in Bppph
  // and restore it at the end of the routine
  // If it is a recipe, we simply construct Vppph in this routine
  auto Bppph( Tcc<TE>::template tensor<Complex<>>("Bppph") );
  auto Vppph( Tcc<TE>::template tensor<Complex<>>("Vppph") );
  auto dressedGph( Tcc<TE>::template tensor<Complex<>>("dressedGph") );
  auto conjGpp( Tcc<TE>::template tensor<Complex<>>("conjGpp") );


  if (arguments->isGiven("coulombVertex")) {
    auto coulombVertex( arguments->getPtr<TensorSet<Complex<>,TE>>("coulombVertex") );
    auto GammaGpp( coulombVertex->get("pp") );
    auto GammaGph( coulombVertex->get("ph") );

    COMPILE(
      (*conjGpp)["Gqr"] <<=
        cc4s::map<Complex<>>(conj<Complex<>>, (*GammaGpp)["Grq"]),
      (*dressedGph)["Gck"] <<= (*GammaGpp)["Gcf"] * (*Tph)["fk"],
      (*Vppph)["bcek"] <<= (*conjGpp)["Gbe"] * (*GammaGph)["Gck"]
    )->execute();
  } else {
    // now coulombIntegrals is a tensor (not a recipe)
    auto whatever(coulombIntegrals->get("ppph"));
    auto recipe(dynamicPtrCast<TensorRecipe<Complex<>,TE>>(whatever));
    if (recipe) {
      std::string msg("No CoulombVertex given but Integrals appear as Recipe. Abort!");
      ASSERT_LOCATION(false, msg, SOURCE_LOCATION);
    }
    COMPILE(
      (*Vppph)["bcek"] <<= (*whatever)["bcek"],
      (*Bppph)["bcek"] <<= (*Vppph)["bcek"]
    )->execute();
  }

  auto pXiajb( Tcc<TE>::template tensor<Complex<>>("pXiajb") );
  auto pXiabj( Tcc<TE>::template tensor<Complex<>>("pXiabj") );
  auto Xijka( Tcc<TE>::template tensor<Complex<>>("Xijka") );
  auto dpXiajb( Tcc<TE>::template tensor<Complex<>>("dpXiajb") );
  auto dpXiabj( Tcc<TE>::template tensor<Complex<>>("dpXiabj") );
  auto Iia( Tcc<TE>::template tensor<Complex<>>("Iia") );
  auto Ippph( Tcc<TE>::template tensor<Complex<>>("Ippph") );
  auto dpImcjk( Tcc<TE>::template tensor<Complex<>>("dpImcjk") );
  auto dpIjkmc( Tcc<TE>::template tensor<Complex<>>("dpIjkmc") );

  std::shared_ptr<TensorSet<Complex<>, TE>> J = std::make_shared<TensorSet<Complex<>, TE>>();
  J->get("ppph") = Vppph;
  J->get("hhhp") = dpIjkmc;



  // Implementing the equations from Piecuch, Computer Physics Communications 149, 71 (2002)
  // Equations in Table 1: we will write the corresponding line number
  // Note that we need only half of the Equations for the intermediates Jppph and Jhhhp
  COMPILE(
    // Line 7: Xijka
    (*Xijka)["ijka"] <<= (*Vhhhp)["ijka"],
    (*Xijka)["ijka"] += (1.0) * (*Tph)["ek"] * (*Vhhpp)["ijea"],
    // Line 6: pXaibc (we store it in a ppph object for memory reasons)
    (*Ippph)["bcai"] <<= (*Vppph)["bcai"]
  )->execute();
  Ippph->sum( Complex<>(1.0), "bcai", "bcai", conj<Complex<>>);

  COMPILE(
    // Line 6: pXaibc
    (*Ippph)["bcai"] += (-0.5) * (*Tph)["am"] * (*Vhhpp)["mibc"],
    // Line 1: Iia
    (*Iia)["ia"] <<= ( 2.0) * (*Vhhpp)["imae"] * (*Tph)["em"],
    (*Iia)["ia"]  += (-1.0) * (*Vhhpp)["imea"] * (*Tph)["em"],
    // (*Iia)["ia"] += Fia["ia"], // TODO: use non-canonical orbitals

    // Line 10: pXiajb
    (*pXiajb)["iajb"] <<= (*Vhphp)["iajb"],
    (*pXiajb)["iajb"] += (-0.5) * (*Vhhhp)["imjb"] * (*Tph)["am"],
    (*pXiajb)["iajb"] += ( 1.0) * (*Tph)["ej"] * (*Ippph)["beai"],

    // Line 14: pXiabj
    (*pXiabj)["iabj"] <<= (*Vhpph)["iabj"],
    (*pXiabj)["iabj"] += (-0.5) * (*Vhhph)["imbj"] * (*Tph)["am"],
    (*pXiabj)["iabj"] += ( 1.0) * (*Ippph)["ebai"] * (*Tph)["ej"],

    // Line 5: Xaibc (store it in the same object to safe memory)
    (*Ippph)["bcai"] += (-0.5) * (*Tph)["am"] * (*Vhhpp)["mibc"],

    // Line 11: ppXiajb
    (*dpXiajb)["iajb"] <<= (*Vhphp)["iajb"],
    (*dpXiajb)["iajb"] += (-1.0) * (*Vhhhp)["imjb"] * (*Tph)["am"],
    (*dpXiajb)["iajb"] += ( 0.5) * (*Tph)["ej"] * (*Ippph)["beai"],

    // Line 15: ppXiabj
    (*dpXiabj)["iabj"] <<= (*Vhpph)["iabj"],
    (*dpXiabj)["iabj"] += (-1.0) * (*Vhhph)["imbj"] * (*Tph)["am"],
    (*dpXiabj)["iabj"] += ( 0.5) * (*Ippph)["ebai"] * (*Tph)["ej"]

  )->execute();

  // Line 17: we use Vppph as scratch for ppIabci
  if (arguments->isGiven("coulombVertex")) {
    COMPILE(
      (*Vppph)["bcek"] += (*conjGpp)["Gbe"] * (*dressedGph)["Gck"]
    )->execute();
  } else {
    COMPILE(
      //(*dpIbcek)["bcek"] += ( 1.0) * (*Vpppp)["bcef"] * (*Tph)["fk"]
      (*Vppph)["bcek"] += ( 1.0) * (*Vpppp)["bcef"] * (*Tph)["fk"]
    )->execute();
  }

  COMPILE(
    // Line 17: 3rd term - 5th term
    (*Vppph)["abci"] += (-1.0) * (*pXiajb)["maic"] * (*Tph)["bm"],
    (*Vppph)["abci"] += (-1.0) * (*Tph)["am"] * (*pXiabj)["mbci"],
    (*Vppph)["abci"] += (-1.0) * (*Iia)["mc"] * (*Tpphh)["abmi"]
  )->execute();
  cc4s::enable_lowmem(true);
  COMPILE(
//    // Line 17: 6th term - 10th term
    (noRing) ? (
      Tcc<TE>::sequence()
    ) : (
      (*Vppph)["abci"] += ( 2.0) * (*Ippph)["ceam"] * (*Tpphh)["beim"],
      Tcc<TE>::sequence()
    ),
    (*Vppph)["abci"] += (-1.0) * (*Ippph)["ceam"] * (*Tpphh)["ebim"],
    (*Vppph)["abci"] += (-1.0) * (*Ippph)["ecam"] * (*Tpphh)["ebmi"],
    (*Vppph)["abci"] += (-1.0) * (*Tpphh)["eaim"] * (*Ippph)["ecbm"]
  )->execute();
  cc4s::enable_lowmem(false);
  COMPILE(
    (*Vppph)["abci"] += ( 1.0) * (*Tpphh)["bamn"] * (*Xijka)["mnic"],
    // Line 19: 1st to 4th term
    (*dpImcjk)["iajk"] <<= (*Vhphh)["iajk"],
    (*dpImcjk)["iajk"] += (-1.0) * (*Vhhhh)["imjk"] * (*Tph)["am"],
    (*dpImcjk)["iajk"] += ( 1.0) * (*dpXiajb)["iaje"] * (*Tph)["ek"],
    (*dpImcjk)["iajk"] += ( 1.0) * (*Tph)["ej"] * (*dpXiabj)["iaek"],
    // Line 19: 5th - 9th term
    (noRing) ? (
      Tcc<TE>::sequence()
    ) : (
      (*dpImcjk)["iajk"] += ( 2.0) * (*Xijka)["imje"] * (*Tpphh)["aekm"],
      Tcc<TE>::sequence()
    ),
    (*dpImcjk)["iajk"] += (-1.0) * (*Xijka)["imje"] * (*Tpphh)["eakm"],
    (*dpImcjk)["iajk"] += (-1.0) * (*Xijka)["mije"] * (*Tpphh)["eamk"],
    (*dpImcjk)["iajk"] += (-1.0) * (*Tpphh)["aemj"] * (*Xijka)["mike"],
    (*dpImcjk)["iajk"] += ( 1.0) * (*Tpphh)["efkj"] * (*Ippph)["efai"],
    // For technical reasons we want the result in the form Vhhhp
    (*dpIjkmc)["jkmc"] <<= map<Complex<>>( conj<Complex<>>, (*dpImcjk)["mcjk"]),
    (*Ippph)["bcek"] <<= (*Vppph)["bcek"]
  )->execute();

  // If Vppph was given as a tensor we have to restore the original result
  if (!arguments->isGiven("coulombVertex")) {
    COMPILE( (*Vppph)["bcek"] <<= (*Bppph)["bcek"] )->execute();
  }


  auto result( New<MapNode>(SOURCE_LOCATION));
  result->setPtr("cTIntermediate", J);
  return result;

}

template <typename F, typename TE>
std::shared_ptr<TensorSet<F, TE>>
cc4s::ct::getPerturbativeCTReference(
  std::shared_ptr<TensorSet<F, TE>> coulombIntegrals,
  std::shared_ptr<TensorSet<F, TE>> amplitudes
) {

  auto Tph ( amplitudes->get("ph") );
  auto Tpphh( amplitudes->get("pphh") );

  auto Vpppp(coulombIntegrals->get("pppp"));
  auto Vphph(coulombIntegrals->get("phph"));
  auto Vhhhh(coulombIntegrals->get("hhhh"));
  auto Vhhhp(coulombIntegrals->get("hhhp"));
  auto Vppph(coulombIntegrals->get("ppph"));
  auto Vhhpp(coulombIntegrals->get("hhpp"));
  auto Vphhh(coulombIntegrals->get("phhh"));
  auto Vhphp(coulombIntegrals->get("hphp"));
  auto Vphhp(coulombIntegrals->get("phhp"));
  auto Vphpp(coulombIntegrals->get("phpp"));
  auto Vhhph(coulombIntegrals->get("hhph"));
  auto Vhpph(coulombIntegrals->get("hpph"));
  auto Vhphh(coulombIntegrals->get("hphh"));

  // Piecuch intermediates
  auto pXiajb( Tcc<TE>::template tensor<F>("pXiajb") );
  auto pXiabj( Tcc<TE>::template tensor<F>("pXiabj") );
  auto pXaibc( Tcc<TE>::template tensor<F>("pXaibc") );
  auto Xaibc( Tcc<TE>::template tensor<F>("Xaibc") );
  auto Xijka( Tcc<TE>::template tensor<F>("Xijka") );
  auto dpXiajb( Tcc<TE>::template tensor<F>("dpXiajb") );
  auto dpXiabj( Tcc<TE>::template tensor<F>("dpXiabj") );
  auto Iia( Tcc<TE>::template tensor<F>("Iia") );

  // Resulting tensors
  std::shared_ptr<TensorSet<F, TE>> result = std::make_shared<TensorSet<F, TE>>();
  auto dpIbcek( Tcc<TE>::template tensor<F>("dpIbcek") );
  auto dpImcjk( Tcc<TE>::template tensor<F>("dpImcjk") );

  auto dpIjkmc( Tcc<TE>::template tensor<F>("dpIjkmc") );
  result->get("ppph") = dpIbcek;
  result->get("hhhp") = dpIjkmc;

  COMPILE(

    // Build Xijka (X-hhhp intermediate)
    (*Xijka)["ijka"] <<= (*Vhhhp)["ijka"],
    (*Xijka)["ijka"] += (1.0) * (*Tph)["fk"] * (*Vhhpp)["ijfa"],

    // Build pXaibc (prime X-phpp intermediate)
    (*pXaibc)["aibc"] <<= (*Vphpp)["aibc"],
    (*pXaibc)["aibc"] += (-0.5) * (*Tph)["al"] * (*Vhhpp)["libc"],

    // Build Iia (I-hp intermediate)
    (*Iia)["ia"] <<= ( 2.0) * (*Vhhpp)["ilaf"] * (*Tph)["fl"],
    (*Iia)["ia"]  += (-1.0) * (*Vhhpp)["ilfa"] * (*Tph)["fl"],
    // (*Iia)["ia"] += Fia["ia"], // TODO: use non-canonical orbitals

    // Build pXiajb (prime X-hphp intermediate)
    (*pXiajb)["iajb"] <<= (*Vhphp)["iajb"],
    (*pXiajb)["iajb"] += (-0.5) * (*Vhhhp)["iljb"] * (*Tph)["al"],
    (*pXiajb)["iajb"] += ( 1.0) * (*Tph)["fj"] * (*pXaibc)["aibf"],

    // Build pXiabj (prime I-hpph intermediate)
    (*pXiabj)["iabj"] <<= (*Vhpph)["iabj"],
    (*pXiabj)["iabj"] += (-0.5) * (*Vhhph)["ilbj"] * (*Tph)["al"],
    (*pXiabj)["iabj"] += ( 1.0) * (*pXaibc)["aifb"] * (*Tph)["fj"],

    // Build Xaibc (X-phpp intermediate)
    (*Xaibc)["aibc"] <<= (*pXaibc)["aibc"],
    (*Xaibc)["aibc"] += (-0.5) * (*Tph)["al"] * (*Vhhpp)["libc"],

    // Build dpXiajb (dp: double prime X-hphp intermediate)
    (*dpXiajb)["iajb"] <<= (*Vhphp)["iajb"],
    (*dpXiajb)["iajb"] += (-1.0) * (*Vhhhp)["iljb"] * (*Tph)["al"],
    (*dpXiajb)["iajb"] += ( 0.5) * (*Tph)["fj"] * (*Xaibc)["aibf"],

    // Build dpXiabj (dp: double prime X-hpph intermediate)
    (*dpXiabj)["iabj"] <<= (*Vhpph)["iabj"],
    (*dpXiabj)["iabj"] += (-1.0) * (*Vhhph)["ilbj"] * (*Tph)["al"],
    (*dpXiabj)["iabj"] += ( 0.5) * (*Xaibc)["aifb"] * (*Tph)["fj"],

    // Build dpIbcek (dp: double prime I-ppph intermediate)
    (*dpIbcek)["bcek"] <<= (*Vppph)["bcek"],
    (*dpIbcek)["bcek"] += ( 1.0) * (*Vpppp)["bcef"] * (*Tph)["fk"],
    (*dpIbcek)["bcek"] += (-1.0) * (*pXiajb)["lbke"] * (*Tph)["cl"],
    (*dpIbcek)["bcek"] += (-1.0) * (*Tph)["bl"] * (*pXiabj)["lcek"],
    (*dpIbcek)["bcek"] += (-1.0) * (*Iia)["le"] * (*Tpphh)["bclk"],
    (*dpIbcek)["bcek"] += ( 2.0) * (*Xaibc)["blef"] * (*Tpphh)["cfkl"],
    (*dpIbcek)["bcek"] += (-1.0) * (*Xaibc)["blef"] * (*Tpphh)["fckl"],
    (*dpIbcek)["bcek"] += (-1.0) * (*Xaibc)["blfe"] * (*Tpphh)["fclk"],
    (*dpIbcek)["bcek"] += (-1.0) * (*Tpphh)["fbkl"] * (*Xaibc)["clfe"],
    (*dpIbcek)["bcek"] += ( 1.0) * (*Tpphh)["cbln"] * (*Xijka)["lnke"],

    // Build dpImcjk (dp: double prime I-hphh intermediate)
    (*dpImcjk)["mcjk"] <<= (*Vhphh)["mcjk"],
    (*dpImcjk)["mcjk"] += (-1.0) * (*Vhhhh)["mljk"] * (*Tph)["cl"],
    (*dpImcjk)["mcjk"] += ( 1.0) * (*dpXiajb)["mcjf"] * (*Tph)["fk"],
    (*dpImcjk)["mcjk"] += ( 1.0) * (*Tph)["fj"] * (*dpXiabj)["mcfk"],
    (*dpImcjk)["mcjk"] += ( 2.0) * (*Xijka)["mljf"] * (*Tpphh)["cfkl"],
    (*dpImcjk)["mcjk"] += (-1.0) * (*Xijka)["mljf"] * (*Tpphh)["fckl"],
    (*dpImcjk)["mcjk"] += (-1.0) * (*Xijka)["lmjf"] * (*Tpphh)["fclk"],
    (*dpImcjk)["mcjk"] += (-1.0) * (*Tpphh)["cflj"] * (*Xijka)["lmkf"],
    (*dpImcjk)["mcjk"] += ( 1.0) * (*Tpphh)["dfkj"] * (*Xaibc)["cmdf"],

    (*dpIjkmc)["jkmc"] <<= map<F>(conj<F>, (*dpImcjk)["mcjk"])
  )->execute();


  return result;
}

template <typename F, typename TE>
double cc4s::ct::getDenominator(
 std::shared_ptr<TensorSet<F, TE>> amplitudes,
 Ptr<cc4s::Tensor<F,TE>> tabcijk,
 Ptr<cc4s::Tensor<F,TE>> zabcijk
) {

  auto Tph ( amplitudes->get("ph") );
  auto Tpphh( amplitudes->get("pphh") );

  // Define the needed quantities for the denominator.
  auto yabcijk( Tcc<TE>::template tensor<F>("yabcijk") );
  auto btabcijk( Tcc<TE>::template tensor<F>("btabcijk") );
  auto Cabij( Tcc<TE>::template tensor<F>("Cabij") );
  auto Wabij( Tcc<TE>::template tensor<F>("Wabij") );
  auto D( Tcc<TE>::template tensor<F>("D") );

  COMPILE(
    // Build y quantity.
    (*yabcijk)["abcijk"] <<= (*Tph)["ai"]
                           * (*Tph)["bj"]
                           * (*Tph)["ck"],
    (*yabcijk)["abcijk"] += (*Tph)["ai"]
                          * (*Tpphh)["bcjk"],
    (*yabcijk)["abcijk"] += (*Tph)["bj"]
                          * (*Tpphh)["acik"],
    (*yabcijk)["abcijk"] += (*Tph)["ck"]
                          * (*Tpphh)["abij"],

    // Build t-bar (T) amplitude.
    (*btabcijk)["abcijk"] <<= (4./3) * (*tabcijk)["abcijk"],
    (*btabcijk)["abcijk"]  += (-2.0) * (*tabcijk)["acbijk"],
    (*btabcijk)["abcijk"]  += (2./3) * (*tabcijk)["bcaijk"],

    // Build C-pphh quantity.
    (*Cabij)["abij"]  <<= (*Tpphh)["abij"],
    (*Cabij)["abij"]   += (*Tph)["ai"] * (*Tph)["bj"],

    // Build W quantity (defined by us) from "antisymmetrized amplitudes".
    (*Wabij)["abij"]  <<=  (2.0) * (*Tpphh)["abij"],
    (*Wabij)["abij"]   += (-1.0) * (*Tpphh)["baij"],

    // Build Piecuch Denominator.
    (*D)[""] <<= (2.0) * (*Tph)["ai"] * map<F>(conj<F>, (*Tph)["ai"]),
    (*D)[""] += (*Wabij)["abij"] * map<F>(conj<F>, (*Cabij)["abij"]),
    (*D)[""] += (*yabcijk)["abcijk"] * map<F>(conj<F>, (*btabcijk)["abcijk"]),
    (*D)[""] += (*yabcijk)["abcijk"] * map<F>(conj<F>, (*zabcijk)["abcijk"])
  )->execute();

  return std::real(D->read()) + 1.0;

}


#define _INSTANTIATE(F, TE) \
  template                  \
  std::shared_ptr< cc4s::TensorSet< F , TE > >   \
  cc4s::ct::getPerturbativeCTReference< F , TE >( \
    std::shared_ptr<cc4s::TensorSet< F , TE > > coulombIntegrals, \
    std::shared_ptr<cc4s::TensorSet< F , TE > > amplitudes \
  );

// Dry tensors
_INSTANTIATE(cc4s::Real<64>, cc4s::DefaultDryTensorEngine)
_INSTANTIATE(cc4s::Complex<64>, cc4s::DefaultDryTensorEngine)
// default tensor engines
_INSTANTIATE(cc4s::Real<64>, cc4s::DefaultTensorEngine)
_INSTANTIATE(cc4s::Complex<64>, cc4s::DefaultTensorEngine)
#undef _INSTANTIATE

#define _INSTANTIATE(F, TE) \
  template                  \
  double   \
  cc4s::ct::getDenominator< F , TE >( \
      std::shared_ptr<TensorSet<F, TE>> amplitudes, \
      Ptr<cc4s::Tensor<F,TE>> tabcijk, \
      Ptr<cc4s::Tensor<F,TE>> zabcijk \
  );
// Dry tensors
_INSTANTIATE(cc4s::Real<64>, cc4s::DefaultDryTensorEngine)
_INSTANTIATE(cc4s::Complex<64>, cc4s::DefaultDryTensorEngine)
// default tensor engines
_INSTANTIATE(cc4s::Real<64>, cc4s::DefaultTensorEngine)
_INSTANTIATE(cc4s::Complex<64>, cc4s::DefaultTensorEngine)

#undef _INSTANTIATE
