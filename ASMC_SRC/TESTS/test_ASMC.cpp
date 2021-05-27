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

#include <utility>
#include <vector>

#include "ASMC.hpp"

#include <Eigen/Core>
#include <fmt/format.h>
#include <fmt/ostream.h>

TEST_CASE("test ASMC decodeAllInJob", "[ASMC]")
{
  DecodingParams params(ASMC_FILE_DIR "/EXAMPLE/exampleFile.n300.array",
                        ASMC_FILE_DIR "/DECODING_QUANTITIES/30-100-2000.decodingQuantities.gz");

  ASMC::ASMC asmc(params);

//  auto result = asmc.decodeAllInJob();
//
//  SECTION("test decode pair summarize")
//  {
//    REQUIRE(result.sumOverPairs.size() == 466440ul);
//  }
}

TEST_CASE("test ASMC decodePairs", "[ASMC]")
{
  ASMC::ASMC asmc(ASMC_FILE_DIR "/EXAMPLE/exampleFile.n300.array",
                  ASMC_FILE_DIR "/DECODING_QUANTITIES/30-100-2000.decodingQuantities.gz");

  std::vector<unsigned> indA = {1, 2, 3};
  std::vector<unsigned> indB = {2, 3, 4};
  auto result = asmc.decodePairs(indA, indB);

  SECTION("test decode pair summarize")
  {
    REQUIRE(result.getPerPairIndices().rows() == 12ll);
    REQUIRE(result.getPerPairIndices().cols() == 4ll);

    REQUIRE(result.getPerPairPosteriorMeans()(11, 0) == Approx(6108.99414f));
    REQUIRE(result.getPerPairPosteriorMeans()(9, 8) == Approx(3151.37231f));
    REQUIRE(result.getPerPairPosteriorMeans()(4, 29) == Approx(27580.66016f));

    REQUIRE(result.getPerPairMAPs()(0, 0) == 29);
    REQUIRE(result.getPerPairMAPs()(9, 1234) == 29);
    REQUIRE(result.getPerPairMAPs()(11, 7) == 22);

    for (Eigen::Index idx = 0ll; idx < result.getPerPairPosteriors().size(); ++idx) {
      REQUIRE((result.getPerPairPosteriors()(idx).colwise().sum() - result.getPerPairPosteriorMeans().row(idx)).isZero(1e-6));
    }
  }
}
