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

#include <test/Test.hpp>
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
  Real<> norm(t.norm2());
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
  Real<> norm(symmetric.norm2());
  REQUIRE(norm == 0.0);
}

