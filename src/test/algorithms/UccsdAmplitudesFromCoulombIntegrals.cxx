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
