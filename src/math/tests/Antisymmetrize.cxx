#include <util/Test.hpp>
#include <ctf.hpp>
#include <Cc4s.hpp>
#include <math/MathFunctions.hpp>

TEST_CASE( "Basic antisymmetrize functionality", "[math]" ) {
  int order(4);
  int lens[] = {2, 2, 2, 2};
  int syms[] = {NS, NS, NS, NS};
  CTF::Tensor<> t(order, lens, syms, CTF::get_universe(), "TestTensor");
  t["aaij"] = 1.0;
  //t.print();
  cc4s::antiSymmetrize("abcd", "bacd", t);
  //t.print();
  double norm(t.norm2());
  REQUIRE(norm == 0.0);
}

TEST_CASE( "Antisymmetrize matrix", "[math]" ) {
  int order(2);
  int lens[] = {2, 2};
  int syms[] = {NS, NS};
  CTF::Tensor<> t(2, lens, syms, CTF::get_universe(), "TestTensor");
  t.fill_random(0.0, 1.0);
  //t.print();
  cc4s::antiSymmetrize("ij", "ji", t);
  //t.print();
  CTF::Scalar<> trace(0.0);
  trace[""] = t["ii"];
  // The trace should be zero
  REQUIRE(trace.get_val() == 0.0);
  // Calculate the symmetric tensor of t It should be zero, and therefore have
  // zero norm
  CTF::Tensor<> symmetric(t);
  symmetric["ij"]  = t["ij"];
  symmetric["ij"] += t["ji"];
  double norm(symmetric.norm2());
  REQUIRE(norm == 0.0);
}

