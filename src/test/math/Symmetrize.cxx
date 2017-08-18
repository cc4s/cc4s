#include <test/Test.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <math/MathFunctions.hpp>

TEST_CASE( "Symmetrize matrix", "[math]" ) {
  int order(2);
  int lens[] = {2, 2};
  int syms[] = {NS, NS};
  CTF::Tensor<> t(2, lens, syms, CTF::get_universe(), "TestTensor");
  t.fill_random(0.0, 1.0);
  t.print();
  cc4s::symmetrize("ij", "ji", t);
  t.print();
  // Calculate the symmetric tensor of t It should be zero, and therefore have
  // zero norm
  t["ij"]  -= t["ji"];
  double norm(t.norm2());
  REQUIRE(norm == 0.0);
}

