#include <test/Test.hpp>

#include <algorithms/UccsdAmplitudesFromCoulombIntegrals.hpp>

#include <util/Log.hpp>

using namespace cc4s;

TEST_CASE( "UccsdAmplitudesFromCoulombIntegrals", "[algorithms]" ) {
  RangeParser range("  2,  3  ,423");
  REQUIRE(range.getRange()[0] == 2);
  REQUIRE(range.getRange()[1] == 3);
  REQUIRE(range.getRange()[2] == 423);
  REQUIRE(range.get_max() == 423);

  RangeParser singleNumber("2");
  REQUIRE(singleNumber.getRange()[0] == 2);
  REQUIRE(singleNumber.get_max() == 2);

  RangeParser singleRange("   0  -  10");
  for (unsigned int i(0) ; i < 11 ; i++) {
    REQUIRE(singleRange.getRange()[i] == i);
  }
  REQUIRE(singleRange.get_max() == 10);

  RangeParser simpleCombination("101,   0  -  10,  2314, 341");
  REQUIRE(simpleCombination.getRange()[0] == 101);
  REQUIRE(simpleCombination.getRange()[12] == 2314);
  REQUIRE(simpleCombination.getRange()[13] == 341);
  for (unsigned int i(1) ; i < 12 ; i++) {
    REQUIRE(simpleCombination.getRange()[i] == i-1);
  }
  REQUIRE(simpleCombination.get_max() == 2314);

}
