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

#include "Data.hpp"
#include "DecodingParams.hpp"
#include "DecodingQuantities.hpp"
#include "FileUtils.hpp"
#include "Individual.hpp"
#include "Types.hpp"
#include <Eigen/Dense>
#include <cstdint>
#include <string>
#include <vector>

using namespace std;

// individual ids and XOR/AND of genotypes
struct PairObservations {
  uint_least8_t iHap;
  uint_least8_t jHap;

  unsigned int iInd;
  unsigned int jInd;

  vector<bool> obsBits;
  vector<bool> homMinorBits;
};

// individual ids and XOR of genotypes
struct DecodingReturnValues {

  /// output for sum over all pairs
  Eigen::ArrayXXf sumOverPairs;

  /// output for sum over all pairs with genotype 00
  Eigen::ArrayXXf sumOverPairs00;

  /// output for sum over all pairs with genotype 01 or 10
  Eigen::ArrayXXf sumOverPairs01;

  /// output for sum over all pairs with genotype 11
  Eigen::ArrayXXf sumOverPairs11;

  int sites = 0;
  unsigned int states = 0;
  vector<bool> siteWasFlippedDuringFolding = {};
};

// does the linear-time decoding
class HMM
{

  int m_batchSize;

  float* m_alphaBuffer;
  float* m_betaBuffer;
  float* m_scalingBuffer;
  float* m_allZeros;

  float* meanPost;
  ushort* MAP;
  float* currentMAPValue;

  // for decoding
  Data& data;
  const DecodingQuantities& m_decodingQuant;
  const DecodingParams decodingParams;

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

  vector<PairObservations> m_observationsBatch;

  // output
  DecodingReturnValues m_decodingReturnValues;
  FileUtils::AutoGzOfstream foutPosteriorMeanPerPair;
  FileUtils::AutoGzOfstream foutMAPPerPair;

public:
  // constructor
  HMM(Data& _data, const DecodingQuantities& _decodingQuant, DecodingParams _decodingParams, bool useBatches,
      int _scalingSkip = 1);

  ~HMM();

  /// Decodes all pairs. Returns a sum of all decoded posteriors (sequenceLength x
  /// states).
  void decodeAll(int jobs, int jobInd);

  vector<vector<float>> decode(const PairObservations& observations);

  pair<vector<float>, vector<float>> decodeSummarize(const PairObservations& observations);

  PairObservations makePairObs(const Individual& iInd, int_least8_t iHap, unsigned int ind1, const Individual& jInd,
                               int_least8_t jHap, unsigned int ind2, int from = 0, int to = 0);

  /// decode a single pair
  ///
  /// i and j must be a valid index in `individuals`
  /// if noBatches is not set then the pair is saved and processing is delayed until the
  /// observationBatch array is full
  ///
  /// @param i index of first individual
  /// @param j index of second individual
  ///
  void decodePair(const uint i, const uint j);

  /// decode a list of pairs
  ///
  /// the input pairs are described using two vectors of indicies of `individuals`.
  /// if noBatches is not set then the pairs are processed in batches for efficiency.
  /// This means that after the call to `decodePairs` there might be unprocessed pairs
  /// waiting for the buffer to be full. Use `finishDecoding` to process these.
  ///
  /// @param individualsA vector of indicies of first individual
  /// @param individualsB vector of indicies of second individual
  ///
  void decodePairs(const vector<uint>& individualsA, const vector<uint>& individualsB);

  /// decode a single pair over a segment of the genome
  ///
  /// i and j must be a valid index in `individuals`, `fromPosition` and `toPosition`
  /// must be less than the size of ?
  /// if noBatches is not set then the pair is saved and processing is delayed until the
  /// observationBatch array is full. This is only efficient if subseqent pairs overlap
  /// in the range `fromPosition` -> `toPosition`
  ///
  /// @param i index of first individual
  /// @param j index of second individual
  ///
  void decodeFromGERMLINE(const uint i, const uint j, const uint fromPosition, const uint toPosition);

  /// returns the current buffer of pair observations
  const vector<PairObservations>& getBatchBuffer()
  {
    return m_observationsBatch;
  }

  /// returns the decoding quantities calculated thus far
  const DecodingReturnValues& getDecodingReturnValues()
  {
    return m_decodingReturnValues;
  }

  /// finish decoding pairs
  ///
  /// tells HMM object to finish processing whatever pairs are stored in the
  /// observationsBatch buffer and close output files
  void finishDecoding();

private:

  void makeBits(PairObservations &obs, unsigned from, unsigned to);

  /// resets the internal state of HMM to a clean state
  void resetDecoding();

  void prepareEmissions();

  // add pair to batch and run if we have enough
  void addToBatch(vector<PairObservations>& obsBatch, const PairObservations& observations);

  // complete with leftover pairs
  void runLastBatch(vector<PairObservations>& obsBatch);

  // decode a batch
  void decodeBatch(const vector<PairObservations>& obsBatch);

  // compute scaling factor for an alpha vector
  void scaleBatch(float* alpha, float* scalings, float* sums, int curBatchSize);

  void applyScaling(float* vec, float* scalings, int curBatchSize);

  // forward step
  void forwardBatch(const float* obsIsZeroBatch, const float* obsIsTwoBatch, int curBatchSize);

  // compute next alpha vector in linear time
  void getNextAlphaBatched(float recDistFromPrevious, float* alphaC, int curBatchSize, const float* previousAlpha,
                           uint pos, const float* obsIsZeroBatch, const float* obsIsTwoBatch, float* AU,
                           float* nextAlpha, const vector<float>& emission1AtSite,
                           const vector<float>& emission0minus1AtSite, const vector<float>& emission2minus0AtSite);
  // backward step
  void backwardBatch(const float* obsIsZeroBatch, const float* obsIsTwoBatch, int curBatchSize);

  // compute previous beta vector in linear time
  void getPreviousBetaBatched(float recDistFromPrevious, int curBatchSize, const float* lastComputedBeta, int pos,
                              const float* obsIsZeroBatch, const float* obsIsTwoBatch, float* vec, float* BU, float* BL,
                              float* currentBeta, const vector<float>& emission1AtSite,
                              const vector<float>& emission0minus1AtSite, const vector<float>& emission2minus0AtSite);

  // --posteriorSums
  void augmentSumOverPairs(vector<PairObservations>& obsBatch, int actualBatchSize, int paddedBatchSize);

  // will eventually write binary output instead of gzipped
  void writePerPairOutput(int actualBatchSize, int paddedBatchSize, const vector<PairObservations>& obsBatch);

  // *********************************************************************
  // non-batched computations (for debugging and pedagogical reasons only)
  // *********************************************************************

  float printVector(const vector<float>& vec);

  float getSumOfVector(const vector<float>& vec);

  vector<float> elementWiseMultVectorScalar(const vector<float>& vec, float val);

  vector<float> elementWiseMultVectorVector(const vector<float>& vec, const vector<float>& factors);

  vector<vector<float>> elementWiseMultMatrixMatrix(const vector<vector<float>>& matrix1,
                                                    const vector<vector<float>>& matrix2);

  vector<vector<float>> normalizeMatrixColumns(const vector<vector<float>>& matrix);

  void fillMatrixColumn(vector<vector<float>>& matrix, const vector<float>& vec, long int pos);

  float roundMorgans(float gen);

  int roundPhysical(int phys);

  vector<float> getEmission(int pos, int distinguished, int undistinguished, int emissionIndex);

  // forward step
  vector<vector<float>> forward(const PairObservations& observations);

  void getNextAlpha(float recDistFromPrevious, vector<float>& alphaC, vector<float>& previousAlpha,
                    vector<float>& nextAlpha, vector<float>& emission1AtSite, vector<float>& emission0minus1AtSite,
                    vector<float>& emission2minus0AtSite, float obsIsZero, float obsIsHomMinor);

  // backward step
  vector<vector<float>> backward(const PairObservations& observations);

  void getPreviousBeta(float recDistFromPrevious, vector<float>& lastComputedBeta, vector<float>& BL, vector<float>& BU,
                       vector<float>& currentBeta, vector<float>& emission1AtSite, vector<float>& emission0minus1AtSite,
                       vector<float>& emission2minus0AtSite, float obsIsZero, float obsIsHomMinor);
};

#endif // ASMC_HMM
