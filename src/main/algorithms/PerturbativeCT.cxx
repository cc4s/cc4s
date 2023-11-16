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
      result = run_tcc<Real<>,TE>(arguments)
    ) || (
      result = run_tcc<Complex<>,TE>(arguments)
    );
  } else {
    using TE = DefaultTensorEngine;
    (
      result = run<Real<>,TE>(arguments)
    ) || (
      result = run<Complex<>,TE>(arguments)
    );
  }
  ASSERT_LOCATION(
    result, "unsupported tensor type as 'coulombIntegrals'",
    arguments->sourceLocation
  );
  return result;
}

template <typename F, typename TE>
Ptr<MapNode> PerturbativeCT::run(
  const Ptr<MapNode> &arguments
) {

  auto amplitudes(
    arguments->getPtr<TensorSet<F,TE>>("amplitudes")
  );
  if (!amplitudes) return nullptr;

  auto tccTph = amplitudes->get("ph");

  auto coulombIntegrals(
    arguments->getPtr<TensorSet<F,TE>>("coulombIntegrals")
  );

  auto
    coulombVertex = arguments->getPtr<TensorSet<Complex<>,TE>>("coulombVertex"),
    GammaGpp = coulombVertex->get("pp");

  auto cTGammaGpp(Tcc<TE>::template tensor<Complex<>>("conjGpp"));
  auto dressedGph(Tcc<TE>::template tensor<Complex<>>("dressedGpp"));
  COMPILE(
    (*cTGammaGpp)["Gqr"] <<=
      cc4s::map<Complex<>>(conj<Complex<>>, (*GammaGpp)["Grq"]),
    (*dressedGph)["Gck"] <<= (*GammaGpp)["Gcf"] * map<Complex<>>(toComplex<F>, (*tccTph)["fk"])
  )->execute();

#define __V__(_idx)                                       \
  ([&arguments]() {                                       \
    COMPILE(arguments                                     \
            ->getPtr<TensorSet<F,TE>>("coulombIntegrals") \
            ->get(#_idx)                                  \
            )->execute();                                 \
    return                                                \
      &(arguments                                         \
        ->getPtr<TensorSet<F,TE>>("coulombIntegrals")     \
        ->get(#_idx)                                      \
        ->evaluate()                                      \
        ->getMachineTensor()                              \
        ->tensor);                                        \
  })()

#define __T__(_idx)                             \
  &(arguments                                   \
    ->getPtr<TensorSet<F, TE>>("amplitudes")    \
    ->get(#_idx)                                \
    ->evaluate()                                \
    ->getMachineTensor()                        \
    ->tensor)

  CTF::Tensor<F>
    &Tph = *__T__(ph),
    &Tpphh = *__T__(phph),
    // Coulomb integrals
    &Vhhhh = *__V__(hhhh),
    &Vhhhp = *__V__(hhhp),
    &Vhhph = *__V__(hhph),
    &Vhhpp = *__V__(hhpp),
    &Vhphh = *__V__(hphh),
    &Vhphp = *__V__(hphp),
    &Vhpph = *__V__(hpph),
    &Vppph = *__V__(ppph);

  int sym[] = {NS, NS, NS, NS};
  int64_t o = Tph.lens[1], v = Tph.lens[0];
  int64_t ooov[] = {o, o, o, v};
  int64_t vovv[] = {v, o, v, v};
  int64_t ovov[] = {o, v, o, v};
  int64_t ovvo[] = {o, v, v, o};
  int64_t ovoo[] = {o, v, o, o};

  CTF::Tensor<F> // Intermediate tensors
    Iia(2, ovov, sym),
    Xijka(4, ooov, sym),
    Xaibc(4, vovv, sym),
    pXiajb(4, ovov, sym),
    pXiabj(4, ovvo, sym),
    dpXiajb(4, ovov, sym),
    dpXiabj(4, ovvo, sym),
    dpImcjk(4, ovoo, sym),
    dpIjkmc(4, ooov, sym);

//TODO
// 1.) Transform T amplitudes, Coulomb Integrals, Coulomb Vertex to CTF
// 2.) we need depending on real/complex different codes for the vertex
// 3.) We have to bring back dpIjkmc and Vppph to tcc objects in the output



  CTF::Univar_Function<F> fConj(conj<F>);

  Xijka["ijka"]  = Vhhhp["ijka"];
  Xijka["ijka"] += Tph["fk"] * Vhhpp["ijfa"];

  Xaibc.sum(1.0, Vppph, "aibc", 0.0, "bcai", fConj);
  Xaibc["aibc"] += (-0.5) * Tph["al"] * Vhhpp["libc"];

  Iia["ia"]  = ( 2.0) * Vhhpp["ilaf"] * Tph["fl"];
  Iia["ia"] += (-1.0) * Vhhpp["ilfa"] * Tph["fl"];

  pXiajb["iajb"]  = Vhphp["iajb"];
  pXiajb["iajb"] += (-0.5) * Vhhhp["iljb"] * Tph["al"];
  pXiajb["iajb"] += ( 1.0) * Tph["fj"] * Xaibc["aibf"];


  pXiabj["iabj"]  = Vhpph["iabj"];
  pXiabj["iabj"] += (-0.5) * Vhhph["ilbj"] * Tph["al"];
  pXiabj["iabj"] += ( 1.0) * Xaibc["aifb"] * Tph["fj"];

  Xaibc["aibc"] += (-0.5) * Tph["al"] * Vhhpp["libc"];

  dpXiajb["iajb"]  = Vhphp["iajb"];
  dpXiajb["iajb"] += (-1.0) * Vhhhp["iljb"] * Tph["al"];
  dpXiajb["iajb"] += ( 0.5) * Tph["fj"] * Xaibc["aibf"];


  dpXiabj["iabj"]  = Vhpph["iabj"];
  dpXiabj["iabj"] += (-1.0) * Vhhph["ilbj"] * Tph["al"];
  dpXiabj["iabj"] += ( 0.5) * Xaibc["aifb"] * Tph["fj"];

  // we use Vppph as scratch for dpIbcek
  //////////////////////////////////////
  CTF::Tensor<Complex<>>
    &_cTGammaGpp = cTGammaGpp->evaluate()->getMachineTensor()->tensor,
    &_dressedGph = dressedGph->evaluate()->getMachineTensor()->tensor;
  if (typeid(Complex<64>) == typeid(F)) {
    Vppph["bcek"] = _cTGammaGpp["Gbe"] * _dressedGph["Gck"];
  } else {

    // ctf transformations for real and imag
    const CTF::Transform<Complex<64>, double>
      real([](Complex<64> d, double &f) { f = std::real(d); }),
      imag([](Complex<64> d, double &f) { f = std::imag(d); });

    int64_t
      gvv[] = {_cTGammaGpp.lens[0], v, v},
      gvo[] = {_cTGammaGpp.lens[0], v, o};

    // these tensors will be once the real part and once the imag part
    CTF::Tensor<F>
      __gamma_Gpp(3, gvv, sym),
      __dressed_Gph(3, gvo, sym);

    real(_cTGammaGpp["Gab"], __gamma_Gpp["Gab"]);
    real(_dressedGph["Gai"], __dressed_Gph["Gai"]);
    Vppph["bcek"]  = __gamma_Gpp["Gbe"] * __dressed_Gph["Gck"];

    imag(_cTGammaGpp["Gab"], __gamma_Gpp["Gab"]);
    imag(_dressedGph["Gai"], __dressed_Gph["Gai"]);
    Vppph["bcek"]  = __gamma_Gpp["Gbe"] * __dressed_Gph["Gck"];

  }

  Vppph["bcek"] += (-1.0) * pXiajb["lbke"] * Tph["cl"];
  Vppph["bcek"] += (-1.0) * Tph["bl"] * pXiabj["lcek"];
  Vppph["bcek"] += (-1.0) * Iia["le"] * Tpphh["bclk"];
  Vppph["bcek"] += ( 2.0) * Xaibc["blef"] * Tpphh["cfkl"];
  Vppph["bcek"] += (-1.0) * Xaibc["blef"] * Tpphh["fckl"];
  Vppph["bcek"] += (-1.0) * Xaibc["blfe"] * Tpphh["fclk"];
  Vppph["bcek"] += (-1.0) * Tpphh["fbkl"] * Xaibc["clfe"];
  Vppph["bcek"] += ( 1.0) * Tpphh["cbln"] * Xijka["lnke"];


  dpImcjk["mcjk"]  = Vhphh["mcjk"];
  dpImcjk["mcjk"] += (-1.0) * Vhhhh["mljk"] * Tph["cl"];
  dpImcjk["mcjk"] += ( 1.0) * dpXiajb["mcjf"] * Tph["fk"];
  dpImcjk["mcjk"] += ( 1.0) * Tph["fj"] * dpXiabj["mcfk"];
  dpImcjk["mcjk"] += ( 2.0) * Xijka["mljf"] * Tpphh["cfkl"];
  dpImcjk["mcjk"] += (-1.0) * Xijka["mljf"] * Tpphh["fckl"];
  dpImcjk["mcjk"] += (-1.0) * Xijka["lmjf"] * Tpphh["fclk"];
  dpImcjk["mcjk"] += (-1.0) * Tpphh["cflj"] * Xijka["lmkf"];
  dpImcjk["mcjk"] += ( 1.0) * Tpphh["dfkj"] * Xaibc["cmdf"];

  dpIjkmc.sum(1.0, dpImcjk, "jkmc", 0.0, "mcjk", fConj);

  using tcc_me = Tensor<F, CtfTensorEngine>;
  auto tcc_tensor = New<tcc_me>(dpImcjk);
  std::shared_ptr<TensorSet<F, TE>> intermediate =
    std::make_shared<TensorSet<F, TE>>();
  intermediate->get("ppph") = New<tcc_me>(Vppph);
  intermediate->get("hhhp") = New<tcc_me>(dpIjkmc);

  auto result( New<MapNode>(SOURCE_LOCATION));
  result->setPtr("cTIntermediate", intermediate);
  return result;

}



template <typename F, typename TE>
Ptr<MapNode> PerturbativeCT::run_tcc(
  const Ptr<MapNode> &arguments
) {

  auto amplitudes(
    arguments->getPtr<TensorSet<F,TE>>("amplitudes")
  );
  if (!amplitudes) return nullptr;

  auto Tph ( amplitudes->get("ph") );
  auto Tpphh( amplitudes->get("pphh") );


  auto coulombIntegrals(
    arguments->getPtr<TensorSet<F,TE>>("coulombIntegrals")
  );

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

  auto coulombVertex(
    arguments->getPtr<TensorSet<Complex<>,TE>>("coulombVertex")
  );
  auto GammaGpp(coulombVertex->get("pp"));

  auto cTGammaGpp(Tcc<TE>::template tensor<Complex<>>("conjGpp"));
  auto dressedGpp(Tcc<TE>::template tensor<Complex<>>("dressedGpp"));
  COMPILE(
    (*cTGammaGpp)["Gqr"] <<=
      cc4s::map<Complex<>>(conj<Complex<>>, (*GammaGpp)["Grq"]),
    (*dressedGpp)["Gck"] <<= (*GammaGpp)["Gcf"] * map<Complex<>>(toComplex<F>, (*Tph)["fk"])
  )->execute();
//    (*Vpppp)["pqsr"] <<= (*cTGammaGpp)["Gps"] * (*GammaGpp)["Gqr"],
//    (*dpIbcek)["bcek"] += ( 1.0) * (*Vpppp)["bcef"] * (*Tph)["fk"],

  // Piecuch intermediates
  auto pXiajb( Tcc<TE>::template tensor<F>("pXiajb") );
  auto pXiabj( Tcc<TE>::template tensor<F>("pXiabj") );
  auto pXaibc( Tcc<TE>::template tensor<F>("pXaibc") );
  auto Xaibc( Tcc<TE>::template tensor<F>("Xaibc") );
  auto Xijka( Tcc<TE>::template tensor<F>("Xijka") );
  auto dpXiajb( Tcc<TE>::template tensor<F>("dpXiajb") );
  auto dpXiabj( Tcc<TE>::template tensor<F>("dpXiabj") );
  auto Iia( Tcc<TE>::template tensor<F>("Iia") );

  auto dpIbcek( Tcc<TE>::template tensor<F>("dpIbcek") );
  auto dpImcjk( Tcc<TE>::template tensor<F>("dpImcjk") );

  auto dpIjkmc( Tcc<TE>::template tensor<F>("dpIjkmc") );

  std::shared_ptr<TensorSet<F, TE>> intermediate =
    std::make_shared<TensorSet<F, TE>>();
  intermediate->get("ppph") = dpIbcek;
  intermediate->get("hhhp") = dpIjkmc;


  COMPILE(

    // Build Xijka (X-hhhp intermediate)
    (*Xijka)["ijka"] <<= (*Vhhhp)["ijka"],
    (*Xijka)["ijka"] += (1.0) * (*Tph)["fk"] * (*Vhhpp)["ijfa"],

    // Build pXaibc (prime X-phpp intermediate)
    (*pXaibc)["aibc"] <<= map<F>(conj<F>, (*Vppph)["bcai"]),
    //(*pXaibc)["aibc"] <<= (*Vphpp)["aibc"],
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

    (*dpIbcek)["bcek"] += map<F>( fromComplex<F>, (*cTGammaGpp)["Gbe"] * (*dressedGpp)["Gck"]),
//    (*dpIbcek)["bcek"] += ( 1.0) * (*Vpppp)["bcef"] * (*Tph)["fk"],
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

    (*dpIjkmc)["jkmc"] <<= map<F>( conj<F>, (*dpImcjk)["mcjk"])
  )->execute();

  auto result( New<MapNode>(SOURCE_LOCATION));
  result->setPtr("cTIntermediate", intermediate);
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
