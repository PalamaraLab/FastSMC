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

#include "HMM.hpp"
#include "FileUtils.hpp"
#include <sstream>

TEST_CASE("test hmm functions", "[HMM]")
{
  DecodingParams params(
      ASMC_FILE_DIR "/EXAMPLE/exampleFile.n300.array",
      ASMC_FILE_DIR "/DECODING_QUANTITIES/30-100-2000.decodingQuantities.gz");
  DecodingQuantities decodingQuantities(params.decodingQuantFile.c_str());
  int sequenceLength = Data::countHapLines(params.hapsFileRoot.c_str());
  Data data(params.hapsFileRoot.c_str(), sequenceLength, decodingQuantities.CSFSSamples,
      params.foldData, params.usingCSFS);
  HMM hmm(data, decodingQuantities, params, !params.noBatches);

  REQUIRE(data.individuals.size() > 20);

  SECTION("test decode pair summarize")
  {
    PairObservations pairObs = makePairObs(data.individuals[0], 1, data.individuals[0], 2);
    vector<vector<float>> decodeResult = hmm.decode(pairObs);
    pair<vector<float>, vector<float>> decodeSummary = hmm.decodeSummarize(pairObs);
    // check that the MAP and posterior mean are the same length
    REQUIRE(decodeSummary.first.size() == decodeSummary.second.size());
    REQUIRE(decodeSummary.first.size() == decodeResult[0].size());
  }

  SECTION("test decode pair")
  {
    REQUIRE(hmm.getBatchBuffer().size() == 0);
    hmm.decodePair(0, 9);
    REQUIRE(hmm.getBatchBuffer().size() == 4);
    for (int i = 0; i < 4; ++i) {
      REQUIRE(hmm.getBatchBuffer()[0].iName == data.individuals[0].name);
      REQUIRE(hmm.getBatchBuffer()[0].jName == data.individuals[9].name);
    }
    hmm.decodePair(1, 1);
    REQUIRE(hmm.getBatchBuffer().size() == 5);
    REQUIRE(hmm.getBatchBuffer()[4].iName == data.individuals[1].name);
    REQUIRE(hmm.getBatchBuffer()[4].jName == data.individuals[1].name);
  }

  SECTION("test decode pairs")
  {
    REQUIRE(hmm.getBatchBuffer().size() == 0);
    hmm.decodePairs({ 0, 1 }, { 9, 1 });
    REQUIRE(hmm.getBatchBuffer().size() == 5);
    for (int i = 0; i < 4; ++i) {
      REQUIRE(hmm.getBatchBuffer()[0].iName == data.individuals[0].name);
      REQUIRE(hmm.getBatchBuffer()[0].jName == data.individuals[9].name);
    }
    REQUIRE(hmm.getBatchBuffer()[4].iName == data.individuals[1].name);
    REQUIRE(hmm.getBatchBuffer()[4].jName == data.individuals[1].name);
  }

  SECTION("test finishDecoding")
  {
    REQUIRE(hmm.getBatchBuffer().size() == 0);
    hmm.decodePair(0, 9);
    REQUIRE(hmm.getBatchBuffer().size() == 4);
    hmm.finishDecoding();
    REQUIRE(hmm.getBatchBuffer().size() == 0);
  }

  SECTION("test fill up buffer")
  {
    // default batch size is 64
    for (int i = 1; i <= 64 / 4; ++i) {
      hmm.decodePair(0, i);
    }

    // buffer should be empty now
    REQUIRE(hmm.getBatchBuffer().size() == 0);
  }
}

TEST_CASE("test hmm with regression test", "[HMM_regression]")
{
  // we only needed to set doPosteriorSums to true, but because C++ does
  // not have keyword arguments we need to go through everything
  DecodingParams params(
      ASMC_FILE_DIR "/EXAMPLE/exampleFile.n300.array",
      ASMC_FILE_DIR "/DECODING_QUANTITIES/30-100-2000.decodingQuantities.gz",
      "", // _outFileRoot
      1, // _jobs
      1, // _jobInd
      "array", // _decodingModeString
      false, // _decodingSequence
      true, // _usingCSFS
      false, // _compress
      false, // _useAncestral
      0.f, // _skipCSFSdistance
      false, // _noBatches
      true // _doPosteriorSums
      );
  DecodingQuantities decodingQuantities(params.decodingQuantFile.c_str());
  int sequenceLength = Data::countHapLines(params.hapsFileRoot.c_str());
  Data data(params.hapsFileRoot.c_str(), sequenceLength, decodingQuantities.CSFSSamples,
      params.foldData, params.usingCSFS);
  HMM hmm(data, decodingQuantities, params, !params.noBatches);

  REQUIRE(data.individuals.size() > 20);

  SECTION("regression test")
  {
    std::string regressionFile = ASMC_FILE_DIR "/../ASMC_SRC/TESTS/data/regression_test_original.gz";
    FileUtils::AutoGzIfstream fin;
    fin.openOrExit(regressionFile);
    hmm.decodeAll(params.jobs, params.jobInd);
    const DecodingReturnValues& decodingReturnValues = hmm.getDecodingReturnValues();
    const vector<vector<float>>& sumOverPairs = decodingReturnValues.sumOverPairs;
    // helpful for debugging, can use these values for the first line,
    // allows us to bypass the decoding step, shows why we need ostringstream
    // vector<vector<float>> sumOverPairs = {{
    //   8.643007,   13.102860,   17.181097,  22.544472,  25.243610,   27.292238,   29.715292,  32.435951,
    //   35.420593,  38.667435,  42.194962,  46.031147,  50.209713,  54.759235,  220.796127,
    //   293.305115,   388.716492, 502.695984, 574.957947, 673.107239,  811.460693, 1020.922913,
    //   1165.540649,   1222.318237,1272.868286, 1159.830933,1104.048706, 1056.524658, 1099.780396,
    //   1469.052979, 788.283264, 811.433533, 838.947144, 865.790283, 891.023132,
    //   916.604370,    943.029663, 966.148865, 985.087158,  993.878296,  990.552917, 976.738892,
    //   958.194580,    940.698242, 925.703125, 912.824768,  900.875122, 888.365906,  874.548523,
    //   859.065186,    842.029297, 823.708496,  804.352112, 784.574463, 764.237488,  743.738098,
    //   723.441650,    704.637756,  686.045715, 665.299805, 642.482788,  617.781250,  591.014709,
    //   561.831421,   529.778381, 493.922882, 452.501678,  400.931580,  311.032745
    // }};
    int pos = 0;
    for( std::string line; getline(fin, line); )
    {
      std::ostringstream oss;
      for (uint k = 0; k < decodingQuantities.states; k++) {
        if (k) oss << "\t";
        oss << sumOverPairs[pos][k];
      }
      REQUIRE(oss.str() == line);
      pos++;
    }
    REQUIRE(pos == sumOverPairs.size());
  }
}
