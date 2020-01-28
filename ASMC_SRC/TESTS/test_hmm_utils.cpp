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

#include <vector>

#include "HmmUtils.hpp"

TEST_CASE("test HMM utility free functions", "[HmmUtils]")
{

  SECTION("test subsetXorVec")
  {
    std::vector<bool> v1 = {false, false, true, true, false, false};
    std::vector<bool> v2 = {false, false, true, true, true, false};

    std::vector<bool> fullXor = {false, false, false, false, true, false};
    REQUIRE(asmc::subsetXorVec(v1, v2) == fullXor);
    REQUIRE(asmc::subsetXorVec(v1, v2, 0, 6) == fullXor);

    std::vector<bool> partialXor = {false, true};
    REQUIRE(asmc::subsetXorVec(v1, v2, 3, 5) == partialXor);
  }

  SECTION("test subsetAndVec")
  {
    std::vector<bool> v1 = {false, false, true, true, false, false};
    std::vector<bool> v2 = {false, false, true, true, true, false};

    std::vector<bool> fullAnd = {false, false, true, true, false, false};
    REQUIRE(asmc::subsetAndVec(v1, v2) == fullAnd);
    REQUIRE(asmc::subsetAndVec(v1, v2, 0, 6) == fullAnd);

    std::vector<bool> partialAnd = {true, false};
    REQUIRE(asmc::subsetAndVec(v1, v2, 3, 5) == partialAnd);
  }
}
