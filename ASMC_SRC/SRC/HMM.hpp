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

#ifndef ASMC_HMM
#define ASMC_HMM

#include "Individual.hpp"
#include "DecodingQuantities.hpp"
#include "DecodingParams.hpp"
#include "Types.hpp"
#include "FileUtils.hpp"
#include "Data.hpp"
#include <string>
#include <vector>

using namespace std;

// individual ids and XOR/AND of genotypes
struct PairObservations {
  string iName, jName;
  int iHap, jHap;
  vector<bool> obsBits;
  vector<bool> homMinorBits;
};

// individual ids and XOR of genotypes
struct DecodingReturnValues {
  vector<vector<float>> sumOverPairs; // output for sum over all pairs
  vector<vector<float>>
      sumOverPairs00; // output for sum over all pairs with genotype 00
  vector<vector<float>>
      sumOverPairs01; // output for sum over all pairs with genotype 01 or 10
  vector<vector<float>>
      sumOverPairs11; // output for sum over all pairs with genotype 11
};

// does the linear-time decoding
class HMM {

  float* alphaBuffer;
  float* betaBuffer;
  float* scalingBuffer;
  float* allZeros;

  // for decoding
  Data& data;
  const DecodingQuantities& decodingQuant;
  const DecodingParams& decodingParams;

  string outFileRoot;
  string expectedCoalTimesFile;

  long int sequenceLength;
  int states;

  int scalingSkip;

  vector<float> expectedCoalTimes;
  vector<bool> useCSFSatThisPosition;

  vector<vector<float>> emission1AtSite;
  vector<vector<float>> emission0minus1AtSite;
  vector<vector<float>> emission2minus0AtSite;

  bool noBatches;
  uint64 currPair = 0;

  const int precision = 2;
  const float minGenetic = 1e-10f;

  // output
  DecodingReturnValues decodingReturnValues;
  FileUtils::AutoGzOfstream foutPosteriorMeanPerPair;
  FileUtils::AutoGzOfstream foutMAPPerPair;
  float* meanPost;
  ushort* MAP;
  float* currentMAPValue;

  public:
  // constructor
  HMM(Data& _data, const DecodingQuantities& _decodingQuant,
      DecodingParams& _decodingParams, bool useBatches, int _scalingSkip = 1);

  void prepareEmissions();

  // Decodes all pairs. Returns a sum of all decoded posteriors (sequenceLength x
  // states).
  DecodingReturnValues decodeAll(int jobs, int jobInd, int batchSize = 64);

  private:
  // add pair to batch and run if we have enough
  void addToBatch(vector<PairObservations>& obsBatch, int batchSize,
      const PairObservations& observations);

  // complete with leftover pairs
  void runLastBatch(vector<PairObservations>& obsBatch);

  // decode a batch
  void decodeBatch(const vector<PairObservations>& obsBatch);

  // compute scaling factor for an alpha vector
  void scaleBatch(float* alpha, float* scalings, float* sums, int curBatchSize);


  void applyScaling(float* vec, float* scalings, int curBatchSize);

  // forward step
  void forwardBatch(
      const float* obsIsZeroBatch, const float* obsIsTwoBatch, int curBatchSize);

  // compute next alpha vector in linear time
  void getNextAlphaBatched(float recDistFromPrevious, float* alphaC, int curBatchSize,
      const float* previousAlpha, uint pos, const float* obsIsZeroBatch,
      const float* obsIsTwoBatch, float* AU, float* nextAlpha,
      const vector<float>& emission1AtSite, const vector<float>& emission0minus1AtSite,
      const vector<float>& emission2minus0AtSite);
  // backward step
  void backwardBatch(
      const float* obsIsZeroBatch, const float* obsIsTwoBatch, int curBatchSize);

  // compute previous beta vector in linear time
  void getPreviousBetaBatched(float recDistFromPrevious, int curBatchSize,
      const float* lastComputedBeta, int pos, const float* obsIsZeroBatch,
      const float* obsIsTwoBatch, float* vec, float* BU, float* BL, float* currentBeta,
      const vector<float>& emission1AtSite, const vector<float>& emission0minus1AtSite,
      const vector<float>& emission2minus0AtSite);

  // --posteriorSums
  void augmentSumOverPairs(
      vector<PairObservations>& obsBatch, int actualBatchSize, int paddedBatchSize);

  // will eventually write binary output instead of gzipped
  void writePerPairOutput(int actualBatchSize, int paddedBatchSize,
      const vector<PairObservations>& obsBatch);

  // *********************************************************************
  // non-batched computations (for debugging and pedagogical reasons only)
  // *********************************************************************

  float printVector(const vector<float>& vec);

  float getSumOfVector(const vector<float>& vec);

  vector<float> elementWiseMultVectorScalar(const vector<float>& vec, float val);

  vector<float> elementWiseMultVectorVector(
      const vector<float>& vec, const vector<float>& factors);

  vector<vector<float>> elementWiseMultMatrixMatrix(
      const vector<vector<float>>& matrix1, const vector<vector<float>>& matrix2);

  vector<vector<float>> normalizeMatrixColumns(const vector<vector<float>>& matrix);

  void fillMatrixColumn(
      vector<vector<float>>& matrix, const vector<float>& vec, long int pos);

  vector<vector<float>> decode(const PairObservations& observations);

  float roundMorgans(float gen);

  int roundPhysical(int phys);

  vector<float> getEmission(
      int pos, int distinguished, int undistinguished, int emissionIndex);

  // forward step
  vector<vector<float>> forward(const PairObservations& observations);

  void getNextAlpha(float recDistFromPrevious, vector<float>& alphaC,
      vector<float>& previousAlpha, vector<float>& nextAlpha,
      vector<float>& emission1AtSite, vector<float>& emission0minus1AtSite,
      vector<float>& emission2minus0AtSite, float obsIsZero, float obsIsHomMinor);

  // backward step
  vector<vector<float>> backward(const PairObservations& observations);

  void getPreviousBeta(float recDistFromPrevious, vector<float>& lastComputedBeta,
      vector<float>& BL, vector<float>& BU, vector<float>& currentBeta,
      vector<float>& emission1AtSite, vector<float>& emission0minus1AtSite,
      vector<float>& emission2minus0AtSite, float obsIsZero, float obsIsHomMinor);
};

#endif // ASMC_HMM
