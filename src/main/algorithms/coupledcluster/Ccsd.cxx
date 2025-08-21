/* Copyright 2021 cc4s.org
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#include <algorithms/coupledcluster/Ccsd.hpp>
#include <MathFunctions.hpp>
#include <Log.hpp>
#include <SharedPointer.hpp>
#include <Exception.hpp>

using namespace cc4s;

template <typename F,typename TE>
CoupledClusterMethodRegistrar<F,TE, Ccsd<F,TE>> Ccsd<F,TE>::registrar_("Ccsd");

template <typename F,typename TE>
std::string Ccsd<F,TE>::describeOptions() {
  std::stringstream stream;
  auto eigenEnergies(
    this->arguments->template getPtr<TensorSet<Real<>,TE>>(
      "slicedEigenEnergies"
    )
  );
  auto epsh(eigenEnergies->get("h"));
  auto No(epsh->inspect()->getLen(0));
  auto integralsSliceSize(
    this->arguments->template getValue<Natural<>>("integralsSliceSize", No)
  );
  stream
    << "integralsSliceSize: " << integralsSliceSize;
  return stream.str();
}

//////////////////////////////////////////////////////////////////////
// Hirata iteration routine for the CCSD amplitudes Tabij and Tai from
// So Hirata, et. al. Chem. Phys. Letters, 345, 475 (2001)
//////////////////////////////////////////////////////////////////////


//TODO: avoid conj in Real version
template <typename F, typename TE>
Ptr<TensorSet<F,TE>> Ccsd<F,TE>::getResiduum(
  const Ptr<TensorSet<F,TE>> &amplitudes
) {
  // construct residuum. Shape will be assumed upon first use.
  auto Rph( Tcc<TE>::template tensor<F>("Rph") );
  auto Rpphh( Tcc<TE>::template tensor<F>("Rpphh") );
  auto residuum(
    New<TensorSet<F,TE>>(
      std::map<std::string,Ptr<TensorExpression<F,TE>>>(
        {{"ph",Rph}, {"pphh",Rpphh}}
      )
    )
  );

  auto coulombIntegrals(
    this->arguments->template getPtr<TensorSet<F,TE>>(
      "coulombIntegrals"
    )
  );
  auto Vpphh(coulombIntegrals->get("pphh"));
  bool ppl(this->arguments->template getValue<bool>("ppl", true));
  bool distinguishable(
    this->arguments->template getValue<int>("distinguishable", 0) == 1
  );

  ASSERT_LOCATION(
    !distinguishable, "distinguishable not supported for complex amplitudes",
    this->arguments->sourceLocation
  );
  if (!amplitudes) {
    // no previous amplitudes given
    COMPILE(
      (*Rph)["ai"] <<= 0.0 * (*Vpphh)["aaii"],
      (*Rpphh)["abij"] <<= (*Vpphh)["abij"]
    )->execute();
  } else {
    // TODO: check if given amplitudes contain expected parts
    // get amplitude parts
    auto Tph( amplitudes->get("ph") );
    auto Tpphh( amplitudes->get("pphh") );
    Tph->inspect()->setName("Tph"); Tpphh->inspect()->setName("Tpphh");

    auto coulombVertex(
      this->arguments->template getPtr<TensorSet<F,TE>>(
        "slicedCoulombVertex"
      )
    );
    auto GammaGpp(coulombVertex->get("pp"));
    auto GammaGph(coulombVertex->get("ph"));
    auto GammaGhp(coulombVertex->get("hp"));
    auto GammaGhh(coulombVertex->get("hh"));

    auto cTGammaGph( Tcc<TE>::template tensor<F>("cTGammaGph"));
    auto cTGammaGhp( Tcc<TE>::template tensor<F>("cTGammaGhp"));
    auto cTGammaGpp( Tcc<TE>::template tensor<F>("cTGammaGpp"));
    auto cTGammaGhh( Tcc<TE>::template tensor<F>("cTGammaGhh"));
    COMPILE(
      (*cTGammaGpp)["Gab"] <<= map<F>(conj<F>, (*GammaGpp)["Gba"]),
      (*cTGammaGhp)["Gia"] <<= map<F>(conj<F>, (*GammaGph)["Gai"]),
      (*cTGammaGph)["Gai"] <<= map<F>(conj<F>, (*GammaGhp)["Gia"]),
      (*cTGammaGhh)["Gij"] <<= map<F>(conj<F>, (*GammaGhh)["Gji"])
    )->execute();
    auto cTDressedGammaGph( Tcc<TE>::template tensor<F>("cTDressedGammaGph"));
    auto cTDressedGammaGpp( Tcc<TE>::template tensor<F>("cTDressedGammaGpp"));
    auto dressedGammaGhh( Tcc<TE>::template tensor<F>("dressedGammaGhh"));
    auto dressedGammaGpp( Tcc<TE>::template tensor<F>("dressedGammaGpp"));

    auto Vphhp(coulombIntegrals->get("phhp"));
    auto Vphph(coulombIntegrals->get("phph"));
    auto Vhhhh(coulombIntegrals->get("hhhh"));
    auto Vhhhp(coulombIntegrals->get("hhhp"));
    auto Vphhh(coulombIntegrals->get("phhh"));
    // Hirata intermediates
    auto Lac( Tcc<TE>::template tensor<F>("Lac") );
    auto Kac( Tcc<TE>::template tensor<F>("Kac") );
    auto Lki( Tcc<TE>::template tensor<F>("Lki") );
    auto Kki( Tcc<TE>::template tensor<F>("Kki") );
    auto Kck( Tcc<TE>::template tensor<F>("Kck") );
    auto Xabij( Tcc<TE>::template tensor<F>("Xabij") ); // T2+T1*T1
    auto Xklij( Tcc<TE>::template tensor<F>("Xklij") );

    // we use Vpphh tensor for Vhhpp. if complex -> conj
    if (std::is_same<F,Complex<>>::value) {
      Vpphh->evaluate()->sum( F(1.0), "abij", "abij", conj<F>);
    }
    {
      auto Xakic( Tcc<TE>::template tensor<F>("Xakic") );

      COMPILE(
        (*Xabij)["abij"] <<= (*Tpphh)["abij"],
        (*Xabij)["abij"] += (*Tph)["ai"] * (*Tph)["bj"],

        (*Kac)["ac"] <<= (-2.0) * (*Vpphh)["cdkl"] * (*Xabij)["adkl"],
        (*Kac)["ac"]  += ( 1.0) * (*Vpphh)["dckl"] * (*Xabij)["adkl"],
        (*Kki)["ki"] <<= ( 2.0) * (*Vpphh)["cdkl"] * (*Xabij)["cdil"],
        (*Kki)["ki"]  += (-1.0) * (*Vpphh)["dckl"] * (*Xabij)["cdil"],

        // Xabij = T2 + 2 * T1*T1
        (*Xabij)["abij"] += (*Tph)["ai"] * (*Tph)["bj"],

        (!distinguishable) ? (
          (*Lac)["ac"] <<= (*Kac)["ac"],
          (*Lki)["ki"] <<= (*Kki)["ki"],
          Tcc<TE>::sequence()
        ) : (
          (*Lac)["ac"] <<= (-1.0) * (*Vpphh)["cdkl"] * (*Xabij)["adkl"],
          (*Lac)["ac"]  += ( 0.5) * (*Vpphh)["dckl"] * (*Xabij)["adkl"],
          (*Lki)["ki"] <<= ( 1.0) * (*Vpphh)["cdkl"] * (*Xabij)["cdil"],
          (*Lki)["ki"]  += (-0.5) * (*Vpphh)["dckl"] * (*Xabij)["cdil"],
          Tcc<TE>::sequence()
        ),
        (*Lac)["ac"] +=  (2.0) * (*cTGammaGpp)["Gac"] * (*GammaGhp)["Gkd"] * (*Tph)["dk"],
        (*Lac)["ac"] += (-1.0) * (*cTGammaGpp)["Gad"] * (*GammaGhp)["Gkc"] * (*Tph)["dk"],

        (*Lki)["ki"] += ( 2.0) * (*Vhhhp)["klic"] * (*Tph)["cl"],
        (*Lki)["ki"] += (-1.0) * (*Vhhhp)["lkic"] * (*Tph)["cl"],

        // Contract Lac with T2 Amplitudes
        (*Rpphh)["abij"] += ( 1.0) * (*Lac)["ac"] * (*Tpphh)["cbij"],

        // Contract Lki with T2 Amplitudes
        (*Rpphh)["abij"] += (-1.0) * (*Lki)["ki"] * (*Tpphh)["abkj"],

        // Contract Coulomb integrals with T1 amplitudes
        (*cTDressedGammaGph)["Gai"] <<= (*cTGammaGph)["Gai"],
        (*cTDressedGammaGph)["Gai"] += (-1.0) * (*cTGammaGhh)["Gki"] * (*Tph)["ak"],
        (*Rpphh)["abij"] += (*cTDressedGammaGph)["Gai"] * (*GammaGpp)["Gbc"] * (*Tph)["cj"],
        (*Rpphh)["abij"] += (-1.0) * (*Vphhh)["akij"] * (*Tph)["bk"],
        (std::is_same<F,Real<>>::value) ? (
          (*Rpphh)["abij"] += (-1.0) * (*Tph)["bk"] * (*Vpphh)["acik"] * (*Tph)["cj"]
        ) : (
          (*Rpphh)["abij"] += (-1.0) * (*Vphhp)["akic"] * (*Tph)["cj"] * (*Tph)["bk"]
        ),

        // Build Xakic
        (*cTDressedGammaGph)["Gai"] += ( 1.0) * (*cTGammaGpp)["Gad"] * (*Tph)["di"],
        (*cTDressedGammaGph)["Gai"] += (-0.5) * (*cTGammaGhp)["Gld"] * (*Xabij)["dail"],
        (*Xakic)["akic"] <<= (*cTDressedGammaGph)["Gai"] * (*GammaGhp)["Gkc"],
        (*Xabij)["dclk"] <<= (1.0) * (*Vpphh)["dclk"],
        (*Xabij)["dclk"] += (-0.5) * (*Vpphh)["cdlk"],
        (*Xakic)["akic"] += (*Xabij)["dclk"] * (*Tpphh)["adil"],
        (*Xabij)["cbkj"] <<= (2.0) * (*Tpphh)["cbkj"],
        (*Xabij)["cbkj"] += (-1.0) * (*Tpphh)["bckj"],
        (*Rpphh)["abij"] += (*Xakic)["akic"] * (*Xabij)["cbkj"]
      )->execute();
    }
    {
      auto Xakci( Tcc<TE>::template tensor<F>("Xakci") );
      COMPILE(
        // Build Xakci
        (*cTDressedGammaGpp)["Gab"] <<= (*cTGammaGpp)["Gab"],
        (*cTDressedGammaGpp)["Gac"] += (-1.0) * (*cTGammaGhp)["Glc"] * (*Tph)["al"],
        (*dressedGammaGhh)["Gij"] <<= (*GammaGhh)["Gij"],
        (*dressedGammaGhh)["Gki"]  += (*GammaGhp)["Gkd"] * (*Tph)["di"],
        // Xakci = Vakci - Vlkci * Tal + Vakcd * Tdi - Vlkcd * Tdail
        (*Xakci)["akci"] <<= (*cTDressedGammaGpp)["Gac"] * (*dressedGammaGhh)["Gki"],
        (*Xakci)["akci"] += (-0.5) * (*Vpphh)["cdlk"] * (*Tpphh)["dail"],
        (*Rpphh)["abij"] += (-1.0) * (*Xakci)["akci"] * (*Tpphh)["cbkj"],
        (*Rpphh)["abij"] += (-1.0) * (*Xakci)["bkci"] * (*Tpphh)["ackj"],

        // Symmetrize Rpphh by applying permutation operator
        (*Rpphh)["abij"] += (*Rpphh)["baji"]

      )->execute();
    }
    //////////////////////////////////////////////////////////////////////
    // Now add all terms to Rpphh that do not need to be symmetrized with
    // the permutation operator
    //////////////////////////////////////////////////////////////////////


    COMPILE(
      // rebuild Xabij = T2 + T1*T1
      (*Xabij)["abij"] <<= (*Tpphh)["abij"],
      (*Xabij)["abij"] += (*Tph)["ai"] * (*Tph)["bj"],

      // Build Xklij intermediate
      (*Xklij)["klij"] <<= (*Vhhhh)["klij"],
      (*Xklij)["klij"]  += (*Vhhhp)["klic"] * (*Tph)["cj"],
      (*Xklij)["klij"]  += (*Vhhhp)["lkjc"] * (*Tph)["ci"],
      // Contract Xklij with T2+T1*T1 Amplitudes via Xabij
      (*Rpphh)["abij"]  += (*Xklij)["klij"] * (*Xabij)["abkl"],
      // Construct last term
      (*Xklij)["klij"] <<= (*Vpphh)["cdkl"] * (*Xabij)["cdij"],
      (ppl) ? (
        (*Rpphh)["abij"] +=  (*Xklij)["klij"] * (*Tpphh)["abkl"]
      ) : (
        (*Rpphh)["abij"] +=  (*Xklij)["klij"] * (*Xabij)["abkl"]
      )
    )->execute();

    if (ppl) {
      Natural<> NG(cTGammaGpp->getLen(0));
      Natural<> Nv(Rpphh->inspect()->getLen(0));
      Natural<> No(Rpphh->inspect()->getLen(2));
      Natural<> sliceSize(
        this->arguments->template getValue<Natural<>>("integralsSliceSize", No)
      );  
      Natural<> numberSlices(Natural<>(ceil(1.0*Nv/sliceSize)));
      std::vector<Ptr<Tensor<F, TE>>> cTSlicedGammaGpp;
      std::vector<Ptr<Tensor<F, TE>>>   SlicedGammaGpp;
      COMPILE(
        (*dressedGammaGpp)["Gab"] <<= (*GammaGpp)["Gab"],
        (*dressedGammaGpp)["Gab"] += (-1.0) * (*GammaGhp)["Gkb"] * (*Tph)["ak"]
      )->execute();
      //Slice GammaGab and store it in a vector
      for (Natural<> v(0); v < numberSlices; v++){
        Natural<> xStart = v*sliceSize;
        Natural<> xEnd = std::min((v+1)*sliceSize,Nv);
        auto dummy(   Tcc<TE>::template tensor<F>("dummy")   );
        auto dummyct( Tcc<TE>::template tensor<F>("dummyct") );
        COMPILE(
          (*dummy )["Gxb"]  <<=
            (*(*dressedGammaGpp)({0, xStart, 0}, {NG, xEnd, Nv}))["Gxb"],
          (*dummyct)["Gxb"] <<=
            (*(*cTDressedGammaGpp)({0, xStart, 0}, {NG, xEnd, Nv}))["Gxb"]
        )->execute();
          SlicedGammaGpp.push_back(dummy);
        cTSlicedGammaGpp.push_back(dummyct);
      }
      // loop over slices
      for (Natural<> m(0); m < numberSlices; m++)
      for (Natural<> n(m); n < numberSlices; n++){
        auto Vxycd( Tcc<TE>::template tensor<F>("Vxycd") );
        auto Rxyij( Tcc<TE>::template tensor<F>("Rxyij") );
        auto Ryxji( Tcc<TE>::template tensor<F>("Ryxji") );
        Natural<> a(n*sliceSize); Natural<> b(m*sliceSize);
        Natural<> Nx(cTSlicedGammaGpp[n]->lens[1]);
        Natural<> Ny(SlicedGammaGpp[m]->lens[1]);
        COMPILE(
          (*Vxycd)["xycd"] <<=
            (*cTSlicedGammaGpp[n])["Gxc"] * (*SlicedGammaGpp[m])["Gyd"],
          (*Rxyij)["xyij"] <<= (*Vxycd)["xycd"] * (*Xabij)["cdij"],
          (*(*Rpphh)({a, b, 0, 0},{a+Nx, b+Ny, No, No}))["xyij"] += (*Rxyij)["xyij"],
          // if a>b: add the same slice at (b,a,j,i):
          (a>b) ? (
            (*Ryxji)["yxji"] <<= (*Rxyij)["xyij"],
            (*(*Rpphh)({b, a, 0, 0},{b+Ny, a+Nx, No, No}))["xyij"] += (*Ryxji)["xyij"]
          ) : (
            Tcc<TE>::sequence()
          )
        )->execute();
      }
    }

    //**************************************************************************
    //***********************  T1 amplitude equations  *************************
    //**************************************************************************
    COMPILE(

      // Contract Kac and Kki with T1 amplitudes
      (*Rph)["ai"] <<= ( 1.0) * (*Kac)["ac"] * (*Tph)["ci"],
      (*Rph)["ai"] +=  (-1.0) * (*Kki)["ki"] * (*Tph)["ak"],

      //Build Kck
      (*Kck)["ck"] <<= ( 2.0) * (*Vpphh)["cdkl"] * (*Tph)["dl"],
      (*Kck)["ck"]  += (-1.0) * (*Vpphh)["dckl"] * (*Tph)["dl"],

      // Contract all the rest terms with T1 and T2 amplitudes
      (*Rph)["ai"] += ( 2.0) * (*Kck)["ck"] * (*Tpphh)["caki"],
      (*Rph)["ai"] += (-1.0) * (*Kck)["ck"] * (*Tpphh)["caik"],
      (*Rph)["ai"] += ( 1.0) * (*Kck)["ck"] * (*Tph)["ci"] * (*Tph)["ak"],
      (std::is_same<F,Real<>>::value) ? (
        (*Rph)["ai"] += ( 2.0) *             (*Vpphh)["acik"]          * (*Tph)["ck"]
      ) : (
        (*Rph)["ai"] += (2.0) * (*cTGammaGph)["Gai"] * (*GammaGhp)["Gkc"] * (*Tph)["ck"]
      ),
      //(*Rph)["ai"] += (-1.0) *           (*Vphph)["akci"]             * (*Tph)["ck"],
      (*Rph)["ai"] += (-1.0) * (*cTGammaGpp)["Gac"] * (*GammaGhh)["Gki"] * (*Tph)["ck"],
      (*Rph)["ai"] +=
         (2.0) * (*cTGammaGpp)["Gac"] * (*GammaGhp)["Gkd"] * (*Xabij)["cdik"],
      (*Rph)["ai"] +=
        (-1.0) * (*cTGammaGpp)["Gad"] * (*GammaGhp)["Gkc"] * (*Xabij)["cdik"],

      (*Rph)["ai"] += (-2.0) * (*Vhhhp)["klic"] * (*Xabij)["ackl"],
      (*Rph)["ai"] += ( 1.0) * (*Vhhhp)["lkic"] * (*Xabij)["ackl"]
    )->execute();
    // If complex: we have to restore Vpphh
    if (std::is_same<F,Complex<>>::value) {
      Vpphh->evaluate()->sum( F(1.0), "abij", "abij", conj<F>);
    }
    COMPILE(
      (*Rpphh)["abij"] += (*Vpphh)["abij"]
    )->execute();
  }
  return residuum;
}

// instantiate
template class cc4s::Ccsd<Real<64>, DefaultDryTensorEngine>;
template class cc4s::Ccsd<Complex<64>, DefaultDryTensorEngine>;
template class cc4s::Ccsd<Real<64>, DefaultTensorEngine>;
template class cc4s::Ccsd<Complex<64>, DefaultTensorEngine>;

