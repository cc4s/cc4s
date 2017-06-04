#include <util/Test.hpp>
#include <math/MathFunctions.hpp>

TEST_CASE( "Test sign of permutation up to degree 4", "[math]" ) {

  // S_2
  REQUIRE(  1 == cc4s::permutationSign("ab", "ab") );
  REQUIRE( -1 == cc4s::permutationSign("ab", "ba") );

  // S_3
  REQUIRE(  1 == cc4s::permutationSign("abc", "abc") );
  REQUIRE( -1 == cc4s::permutationSign("abc", "bac") );
  REQUIRE( -1 == cc4s::permutationSign("abc", "acb") );
  REQUIRE( -1 == cc4s::permutationSign("abc", "cba") );
  REQUIRE(  1 == cc4s::permutationSign("abc", "bca") );
  REQUIRE(  1 == cc4s::permutationSign("abc", "cab") );

  // S_4
  REQUIRE(  1 == cc4s::permutationSign("abcd", "abcd") );
  REQUIRE( -1 == cc4s::permutationSign("abcd", "abdc") );
  REQUIRE( -1 == cc4s::permutationSign("abcd", "acbd") );
  REQUIRE(  1 == cc4s::permutationSign("abcd", "acdb") );
  REQUIRE(  1 == cc4s::permutationSign("abcd", "adbc") );
  REQUIRE( -1 == cc4s::permutationSign("abcd", "adcb") );
  REQUIRE( -1 == cc4s::permutationSign("abcd", "bacd") );
  REQUIRE(  1 == cc4s::permutationSign("abcd", "badc") );
  REQUIRE(  1 == cc4s::permutationSign("abcd", "bcad") );
  REQUIRE( -1 == cc4s::permutationSign("abcd", "bcda") );
  REQUIRE( -1 == cc4s::permutationSign("abcd", "bdac") );
  REQUIRE(  1 == cc4s::permutationSign("abcd", "bdca") );
  REQUIRE(  1 == cc4s::permutationSign("abcd", "cabd") );
  REQUIRE( -1 == cc4s::permutationSign("abcd", "cadb") );
  REQUIRE( -1 == cc4s::permutationSign("abcd", "cbad") );
  REQUIRE(  1 == cc4s::permutationSign("abcd", "cbda") );
  REQUIRE(  1 == cc4s::permutationSign("abcd", "cdab") );
  REQUIRE( -1 == cc4s::permutationSign("abcd", "cdba") );
  REQUIRE( -1 == cc4s::permutationSign("abcd", "dabc") );
  REQUIRE(  1 == cc4s::permutationSign("abcd", "dacb") );
  REQUIRE(  1 == cc4s::permutationSign("abcd", "dbac") );
  REQUIRE( -1 == cc4s::permutationSign("abcd", "dbca") );
  REQUIRE( -1 == cc4s::permutationSign("abcd", "dcab") );
  REQUIRE(  1 == cc4s::permutationSign("abcd", "dcba") );

  // S_4 repeating
  REQUIRE(   1 == cc4s::permutationSign("aacd", "aacd") );
  REQUIRE(  -1 == cc4s::permutationSign("aacd", "aadc") );

}
