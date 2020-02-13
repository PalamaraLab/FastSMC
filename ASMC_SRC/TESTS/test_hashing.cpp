//    This file is part of ASMC, developed by Pier Francesco Palamara.
//
//    ASMC is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    ASMC is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with ASMC.  If not, see <https://www.gnu.org/licenses/>.

#include "catch.hpp"

#include "HASHING/Individuals.hpp"

TEST_CASE("test individuals", "[HASHING]")
{
  Individuals<8, 3> ind(5u);

  REQUIRE(ind.getIdNum() == 5ul);

  // Check up to 10 - but internally we're just going 0-1-2-0-1-2-0-1-2-0
  for (auto i = 0; i < 10; ++i) {
    REQUIRE(ind.getWordHash(i) == 0ul);
    REQUIRE(ind.getWordString(i) == "00000000");
  }

  ind.setMarker(0, 0);
  ind.setMarker(1, 2);

  ind.setMarker(2, 2);
  ind.setMarker(2, 3);

  REQUIRE(ind.getWordHash(0) == 1ul);
  REQUIRE(ind.getWordString(0) == "00000001");

  REQUIRE(ind.getWordHash(1) == 4ul);
  REQUIRE(ind.getWordString(1) == "00000100");

  REQUIRE(ind.getWordHash(2) == 12ul);
  REQUIRE(ind.getWordString(2) == "00001100");

  // Clear 2
  ind.clear(2);
  REQUIRE(ind.getWordHash(2) == 0ul);
  REQUIRE(ind.getWordString(2) == "00000000");

  // Clear 1 by clearing 4
  ind.clear(4);
  REQUIRE(ind.getWordHash(1) == 0ul);
  REQUIRE(ind.getWordString(1) == "00000000");
}
