#include <test/Test.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <math/FockVector.hpp>

TEST_CASE( "Basic FockVector testing", "[math]" ) {

  int No = 10;
  int Nv = 4;
  int vo[] = {Nv, No};
  int vvoo[] = {Nv, Nv, No, No};
  int ov[] = {No, Nv};
  int oovv[] = {No, No, Nv, Nv};

  int symsSingles[] = {NS, NS};
  int symsDoubles[] = {NS, NS, NS, NS};

  CTF::Tensor<double> Rai(2, vo, symsSingles, CTF::get_universe(), "Rai");
  CTF::Tensor<double> Rabij(
    4, vvoo, symsDoubles, CTF::get_universe(), "Rabij");
  cc4s::FockVector<double> R(
    {Rai, Rabij}, {"ai", "abij"}
  );

  CTF::Tensor<double> Lia(2, vo, symsSingles, CTF::get_universe(), "Lia");
  CTF::Tensor<double> Lijab(
    4, vvoo, symsDoubles, CTF::get_universe(), "Lijab");
  cc4s::FockVector<double> L(
    {Lia, Lijab}, {"ia", "ijab"}
  );

  double lNorm(L.dot(L));
  REQUIRE(lNorm == 0.0);

  double rNorm(R.dot(R));
  REQUIRE(rNorm == 0.0);

}

