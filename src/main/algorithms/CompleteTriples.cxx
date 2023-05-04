#include <memory>

#include <algorithms/CompleteTriples.hpp>
#include <TensorSet.hpp>

using namespace cc4s;


template <typename F, typename TE>
std::shared_ptr<TensorSet<F, TE>>
cc4s::ct::getCompleteTriples(
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
  auto Vpphp(coulombIntegrals->get("pphp"));
  auto Vphhh(coulombIntegrals->get("phhh"));
  auto Vhphp(coulombIntegrals->get("hphp"));
  auto Vphhp(coulombIntegrals->get("phhp"));
  auto Vphpp(coulombIntegrals->get("phpp"));
  auto Vhhph(coulombIntegrals->get("hhph"));
  auto Vhppp(coulombIntegrals->get("hppp"));
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
  result->get("ppph") = dpIbcek;
  result->get("hphh") = dpImcjk;

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
    (*dpImcjk)["mcjk"] += ( 1.0) * (*Tpphh)["dfkj"] * (*Xaibc)["cmdf"]

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
  (*D)[""] <<= (2.0) * (*Tph)["ai"] * (*Tph)["ai"],
  (*D)[""] += (*Wabij)["abij"] * (*Cabij)["abij"],
  (*D)[""] += (*yabcijk)["abcijk"] * (*btabcijk)["abcijk"],
  (*D)[""] += (*yabcijk)["abcijk"] * (*zabcijk)["abcijk"]
  
  )->execute();
  
  return std::real(D->read()) + 1.0;

}


#define _INSTANTIATE(F, TE) \
  template                  \
  std::shared_ptr< cc4s::TensorSet< F , TE > >   \
  cc4s::ct::getCompleteTriples< F , TE >( \
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
