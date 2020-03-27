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

#include <iostream>
#include <string>
#include <vector>

#include "Data.hpp"
#include "DecodingParams.hpp"
#include "DecodingQuantities.hpp"
#include "FastSMC.hpp"
#include "FileUtils.hpp"
#include "HMM.hpp"
#include "Timer.hpp"



TEST_CASE("test FastSMC HMM with regression test", "[FastSMC_regression]")
{
  srand(1234);

  DecodingParams params;
  params.decodingQuantFile = ASMC_FILE_DIR
      "/FASTSMC_EXAMPLE/out.25.n300.chr2.len30.dens1.disc10-20-2000.demoCEU.mapnorm.array.decodingQuantities.gz";
  params.inFileRoot =
      ASMC_FILE_DIR "/FASTSMC_EXAMPLE/out.25.n300.chr2.len30.dens1.disc10-20-2000.demoCEU.mapnorm.array";
  params.outFileRoot = "/tmp/FastSMCresults";
  params.decodingModeString = "array";
  params.decodingMode = DecodingMode::arrayFolded;
  params.foldData = true;
  params.usingCSFS = true;
  params.batchSize = 32;
  params.recallThreshold = 3;
  params.min_m = 1.5;
  params.GERMLINE = true;
  params.BIN_OUT = false;
  params.time = 50;
  params.noConditionalAgeEstimates = true;
  params.doPerPairMAP = true;
  params.doPerPairPosteriorMean = true;

  // read decoding quantities from file
  DecodingQuantities decodingQuantities(params.decodingQuantFile.c_str());
  const auto sequenceLength = Data::countHapLines(params.inFileRoot);

  Data data(params.inFileRoot, sequenceLength, decodingQuantities.CSFSSamples, params.foldData, params.usingCSFS,
            params.jobInd, params.jobs);

  HMM hmm(data, decodingQuantities, params, !params.noBatches);
  hmm.decodeAll(params.jobs, params.jobInd);

  ASMC::RunFastSMC(params, data, hmm);

  SECTION("regression test")
  {
    const auto expectedNumLines = 1512ul;

    // Read lines from existing regression test output into a vector of strings
    std::vector<std::string> regressionLines;
    regressionLines.reserve(expectedNumLines);
    {
      FileUtils::AutoGzIfstream fin_regression;
      fin_regression.openOrExit(ASMC_FILE_DIR "/FASTSMC_EXAMPLE/regression_output.ibd.gz");
      for (std::string f_line; getline(fin_regression, f_line);) {
        regressionLines.emplace_back(f_line);
      }
      fin_regression.close();
    }

    // Read lines from generated test output into a vector of strings
    std::vector<std::string> generatedLines;
    generatedLines.reserve(expectedNumLines);
    {
      FileUtils::AutoGzIfstream fin_generated;
      fin_generated.openOrExit(params.outFileRoot + ".1.1.FastSMC.ibd.gz");
      for (std::string f_line; getline(fin_generated, f_line);) {
        generatedLines.emplace_back(f_line);
      }
      fin_generated.close();
    }

    REQUIRE(regressionLines.size() == expectedNumLines);
    REQUIRE(regressionLines.size() == generatedLines.size());

    for (auto lineNum = 0ul; lineNum < regressionLines.size(); ++lineNum) {
      REQUIRE(regressionLines.at(lineNum) == generatedLines.at(lineNum));
    }
  }
}
