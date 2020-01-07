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

#include "HMM.hpp"

#include <cassert>
#include <exception>
#include <iostream>
#include <limits>
#include <map>

#include <emmintrin.h>
#include <math.h>
#include <pmmintrin.h>
#include <xmmintrin.h>

#include "MemoryUtils.hpp"
#include "Timer.hpp"
//#include <sys/types.h>

#include "StringUtils.hpp"

#include <Eigen/Dense>
#include <chrono>
#include <sstream>
#include <utility>

using namespace std;

#ifdef NO_SSE
#define MODE "NO_SSE"
#define VECX 4
#endif

// SSE vectorization (block size = 4)
#ifdef SSE
#define MODE "SSE"
#define VECX 4
#define FLOAT __m128
#define LOAD _mm_load_ps
#define STORE _mm_store_ps
#define MULT _mm_mul_ps
#define ADD _mm_add_ps
#define RECIPROCAL _mm_rcp_ps
#define LOAD1 _mm_load1_ps
#endif

// AVX vectorization (block size = 8)
#ifdef AVX
#define MODE "AVX"
#include <immintrin.h>
#define VECX 8
#define FLOAT __m256
#define LOAD _mm256_load_ps
#define STORE _mm256_store_ps
#define MULT _mm256_mul_ps
#define ADD _mm256_add_ps
#define RECIPROCAL _mm256_rcp_ps
#define LOAD1 _mm256_broadcast_ss
#endif

// AVX512 vectorization (block size = 16)
#ifdef AVX512
#define MODE "AVX512"
#include <immintrin.h>
#define VECX 16
#define FLOAT __m512
#define LOAD _mm512_load_ps
#define STORE _mm512_store_ps
#define MULT _mm512_mul_ps
#define ADD _mm512_add_ps
#define RECIPROCAL _mm512_rcp14_ps
#define LOAD1 _mm512_set1_ps
#endif

// to keep track of time
std::chrono::high_resolution_clock::duration t1sum(0), t2sum(0), t1sumBack(0);
std::chrono::high_resolution_clock::duration ticksForward(0), ticksBackward(0), ticksCombine(0), ticksSumOverPairs(0),
    ticksOutputPerPair(0);

// gets genotypes for decoding  (xor --> 1 if het)
vector<bool> xorVec(const vector<bool>& x, const vector<bool>& y, const unsigned from = 0,
                    const unsigned to = numeric_limits<unsigned>::max()) noexcept
{
  const unsigned min_to = std::min(static_cast<unsigned>(x.size()), to);

  assert(x.size() == y.size());
  assert(from < min_to);

  vector<bool> ret(min_to - from);

  for (unsigned i = from; i < min_to; i++) {
    ret[i - from] = x[i] ^ y[i];
  }

  return ret;
}

// computes and of genotype, used to distinguish homozygous minor/derived from
// homozygous major/ancestral
vector<bool> andVec(const vector<bool>& x, const vector<bool>& y, const unsigned from = 0,
                    const unsigned to = numeric_limits<unsigned>::max()) noexcept
{
  const unsigned min_to = std::min(static_cast<unsigned>(x.size()), to);

  assert(x.size() == y.size());
  assert(from < min_to);

  vector<bool> ret(min_to - from);

  for (unsigned i = from; i < min_to; i++) {
    ret[i - from] = x[i] & y[i];
  }

  return ret;
}

// read expected times from a file
vector<float> readExpectedTimesFromIntervalsFil(const char* fileName)
{
  FileUtils::AutoGzIfstream br;
  br.openOrExit(fileName);
  string line;
  vector<float> expCoalTimes;
  while (getline(br, line)) {
    vector<string> splitString;
    istringstream iss(line);
    string buf;
    while (iss >> buf)
      splitString.push_back(buf);
    if (splitString.size() != 3) {
      cerr << fileName
           << " should have \"intervalStart\texpectedCoalescentTime\tintervalEnd\" at "
              "each line."
           << endl;
      exit(1);
    }
    expCoalTimes.push_back(StringUtils::stof(splitString[1]));
  }
  return expCoalTimes;
}

// to print time
void printPctTime(const char* str, double fracTime)
{
  printf("Time in %-15s: %4.1f%%\n", str, 100 * fracTime);
}

// constructor
HMM::HMM(Data& _data, const DecodingQuantities& _decodingQuant, DecodingParams _decodingParams, bool useBatches,
         int _scalingSkip)
    : m_batchSize(64), data(_data), m_decodingQuant(_decodingQuant), decodingParams(_decodingParams),
      scalingSkip(_scalingSkip), noBatches(!useBatches)
{

  cout << "Will decode using " << MODE << " instruction set.\n\n";
  outFileRoot = decodingParams.outFileRoot;
  expectedCoalTimesFile = decodingParams.expectedCoalTimesFile;
  sequenceLength = data.sites;
  states = _decodingQuant.states;
  useCSFSatThisPosition = vector<bool>(sequenceLength, false);
  emission1AtSite = vector<vector<float>>(sequenceLength, vector<float>(states));
  emission0minus1AtSite = vector<vector<float>>(sequenceLength, vector<float>(states));
  emission2minus0AtSite = vector<vector<float>>(sequenceLength, vector<float>(states));
  prepareEmissions();

  for (int i = 0; i < m_batchSize; i++) {
    fromBatch.push_back(0);
    toBatch.push_back(sequenceLength);
  }

  // get state threshold
  stateThreshold = getStateThreshold();
  // probabilityThreshold = (1./decodingQuant->states)*stateThreshold;
  probabilityThreshold = 0;
  for (int i = 0; i < stateThreshold; i++) {
    probabilityThreshold += m_decodingQuant.initialStateProb[i];
  }

  if (decodingParams.noConditionalAgeEstimates) {
    ageThreshold = states;
  } else {
    ageThreshold = stateThreshold;
  }

  startBatch = sequenceLength;
  endBatch = 0;

  if (decodingParams.doPerPairPosteriorMean) {
    expectedCoalTimes = readExpectedTimesFromIntervalsFil(expectedCoalTimesFile.c_str());
  }

  // allocate buffers
  m_scalingBuffer = ALIGNED_MALLOC_FLOATS(m_batchSize);
  m_alphaBuffer = ALIGNED_MALLOC_FLOATS(sequenceLength * states * m_batchSize);
  m_betaBuffer = ALIGNED_MALLOC_FLOATS(sequenceLength * states * m_batchSize);
  m_allZeros = ALIGNED_MALLOC_FLOATS(sequenceLength * m_batchSize);
  memset(m_allZeros, 0, sequenceLength * m_batchSize * sizeof(m_allZeros[0]));

  if (decodingParams.doPerPairPosteriorMean) {
    meanPost = ALIGNED_MALLOC_FLOATS(sequenceLength * m_batchSize);
  }
  if (decodingParams.doPerPairMAP) {
    MAP = ALIGNED_MALLOC_USHORTS(sequenceLength * m_batchSize);
    currentMAPValue = ALIGNED_MALLOC_FLOATS(m_batchSize);
  }

  resetDecoding();

  // output for python interface (TODO: not sure if this is the right place)
  m_decodingReturnValues.sites = _data.sites;
  m_decodingReturnValues.states = _decodingQuant.states;
  m_decodingReturnValues.siteWasFlippedDuringFolding = _data.siteWasFlippedDuringFolding;
}

HMM::~HMM()
{
  ALIGNED_FREE(m_betaBuffer);
  ALIGNED_FREE(m_alphaBuffer);
  ALIGNED_FREE(m_scalingBuffer);
  ALIGNED_FREE(m_allZeros);

  if (decodingParams.doPerPairPosteriorMean) {
    ALIGNED_FREE(meanPost);
  }
  if (decodingParams.doPerPairMAP) {
    ALIGNED_FREE(MAP);
    ALIGNED_FREE(currentMAPValue);
  }
}

PairObservations HMM::makePairObs(const Individual& iInd, int_least8_t iHap, unsigned int ind1, const Individual& jInd,
                                  int_least8_t jHap, unsigned int ind2, int from, int to)
{
  PairObservations ret;
  ret.iHap = iHap;
  ret.jHap = jHap;
  ret.iInd = ind1;
  ret.jInd = ind2;
  if (!decodingParams.GERMLINE || noBatches) {
    ret.obsBits = xorVec(iHap == 1 ? iInd.genotype1 : iInd.genotype2, jHap == 1 ? jInd.genotype1 : jInd.genotype2, 0,
                         sequenceLength);
    ret.homMinorBits = andVec(iHap == 1 ? iInd.genotype1 : iInd.genotype2, jHap == 1 ? jInd.genotype1 : jInd.genotype2,
                              0, sequenceLength);
  }
  return ret;
}

void HMM::makeBits(PairObservations& obs, unsigned from, unsigned to)
{
  unsigned iInd = obs.iInd;
  unsigned jInd = obs.jInd;
  obs.obsBits = xorVec(obs.iHap == 1 ? data.individuals[iInd].genotype1 : data.individuals[iInd].genotype2,
                       obs.jHap == 1 ? data.individuals[jInd].genotype1 : data.individuals[jInd].genotype2, from, to);
  obs.homMinorBits =
      andVec(obs.iHap == 1 ? data.individuals[iInd].genotype1 : data.individuals[iInd].genotype2,
             obs.jHap == 1 ? data.individuals[jInd].genotype1 : data.individuals[jInd].genotype2, from, to);
}

void HMM::prepareEmissions()
{
  if (decodingParams.skipCSFSdistance < std::numeric_limits<float>::infinity()) {
    useCSFSatThisPosition[0] = true;
    float lastGenCSFSwasUsed = 0.f;
    for (int pos = 1; pos < sequenceLength; pos++) {
      if (data.geneticPositions[pos] - lastGenCSFSwasUsed >= decodingParams.skipCSFSdistance) {
        // this position is a CSFS position
        useCSFSatThisPosition[pos] = true;
        lastGenCSFSwasUsed = data.geneticPositions[pos];
      }
    }
  }
  for (int pos = 0; pos < sequenceLength; pos++) {
    if (useCSFSatThisPosition[pos]) {
      int undistAtThisSiteFor0dist = data.undistinguishedCounts[pos][0];
      int undistAtThisSiteFor1dist = data.undistinguishedCounts[pos][1];
      int undistAtThisSiteFor2dist = data.undistinguishedCounts[pos][2];
      if (decodingParams.foldData) {
        // working with folded data
        for (int k = 0; k < states; k++) {
          if (undistAtThisSiteFor1dist >= 0) {
            emission1AtSite[pos][k] = decodingParams.decodingSequence
                                          ? m_decodingQuant.foldedCSFSmap[undistAtThisSiteFor1dist][1][k]
                                          : m_decodingQuant.foldedAscertainedCSFSmap[undistAtThisSiteFor1dist][1][k];
          } else {
            emission1AtSite[pos][k] = 0.f;
          }
          emission0minus1AtSite[pos][k] =
              decodingParams.decodingSequence
                  ? (m_decodingQuant.foldedCSFSmap[undistAtThisSiteFor0dist][0][k] - emission1AtSite[pos][k])
                  : (m_decodingQuant.foldedAscertainedCSFSmap[undistAtThisSiteFor0dist][0][k] -
                     emission1AtSite[pos][k]);
          if (undistAtThisSiteFor2dist >= 0) {
            emission2minus0AtSite[pos][k] =
                decodingParams.decodingSequence
                    ? (m_decodingQuant.foldedCSFSmap[undistAtThisSiteFor2dist][0][k] -
                       m_decodingQuant.foldedCSFSmap[undistAtThisSiteFor0dist][0][k])
                    : (m_decodingQuant.foldedAscertainedCSFSmap[undistAtThisSiteFor2dist][0][k] -
                       m_decodingQuant.foldedAscertainedCSFSmap[undistAtThisSiteFor0dist][0][k]);
          } else {
            emission2minus0AtSite[pos][k] =
                decodingParams.decodingSequence
                    ? (0 - m_decodingQuant.foldedCSFSmap[undistAtThisSiteFor0dist][0][k])
                    : (0 - m_decodingQuant.foldedAscertainedCSFSmap[undistAtThisSiteFor0dist][0][k]);
          }
        }
      } else {
        // working with unfolded data
        for (int k = 0; k < states; k++) {
          if (undistAtThisSiteFor1dist >= 0) {
            emission1AtSite[pos][k] = decodingParams.decodingSequence
                                          ? m_decodingQuant.CSFSmap[undistAtThisSiteFor1dist][1][k]
                                          : m_decodingQuant.ascertainedCSFSmap[undistAtThisSiteFor1dist][1][k];
          } else {
            emission1AtSite[pos][k] = 0.f;
          }
          float emission0AtThisSiteAndState = 0.f;
          if (undistAtThisSiteFor0dist >= 0) {
            emission0AtThisSiteAndState = decodingParams.decodingSequence
                                              ? m_decodingQuant.CSFSmap[undistAtThisSiteFor0dist][0][k]
                                              : m_decodingQuant.ascertainedCSFSmap[undistAtThisSiteFor0dist][0][k];
          }
          emission0minus1AtSite[pos][k] = emission0AtThisSiteAndState - emission1AtSite[pos][k];
          if (undistAtThisSiteFor2dist >= 0) {
            int dist = 2;
            int undist = undistAtThisSiteFor2dist;
            if (undistAtThisSiteFor2dist == data.totalSamplesBound - 2) {
              // for monomorphic derived, fold to CSFS[0][0]
              dist = 0;
              undist = 0;
            }
            emission2minus0AtSite[pos][k] =
                decodingParams.decodingSequence
                    ? (m_decodingQuant.CSFSmap[undist][dist][k] - emission0AtThisSiteAndState)
                    : (m_decodingQuant.ascertainedCSFSmap[undist][dist][k] - emission0AtThisSiteAndState);
          } else {
            emission2minus0AtSite[pos][k] = 0 - emission0AtThisSiteAndState;
          }
        }
      }
    } else {
      // this position is not a CSFS position
      for (int k = 0; k < states; k++) {
        emission1AtSite[pos][k] = decodingParams.decodingSequence ? m_decodingQuant.classicEmissionTable[1][k]
                                                                  : m_decodingQuant.compressedEmissionTable[1][k];
        emission0minus1AtSite[pos][k] =
            decodingParams.decodingSequence
                ? (m_decodingQuant.classicEmissionTable[0][k] - m_decodingQuant.classicEmissionTable[1][k])
                : (m_decodingQuant.compressedEmissionTable[0][k] - m_decodingQuant.compressedEmissionTable[1][k]);
        // emission2 = emission0
        emission2minus0AtSite[pos][k] = 0.f;
      }
    }
  }
}

void HMM::resetDecoding()
{
  if (decodingParams.doPerPairPosteriorMean) {
    if (foutPosteriorMeanPerPair) {
      foutPosteriorMeanPerPair.close();
    }
    foutPosteriorMeanPerPair.openOrExit(outFileRoot + ".perPairPosteriorMeans.gz");
  }
  if (decodingParams.doPerPairMAP) {
    if (foutMAPPerPair) {
      foutMAPPerPair.close();
    }
    foutMAPPerPair.openOrExit(outFileRoot + ".perPairMAP.gz");
  }

  m_decodingReturnValues.sumOverPairs = Eigen::ArrayXXf::Zero(sequenceLength, states);

  if (decodingParams.doMajorMinorPosteriorSums) {
    m_decodingReturnValues.sumOverPairs00 = Eigen::ArrayXXf::Zero(sequenceLength, states);
    m_decodingReturnValues.sumOverPairs01 = Eigen::ArrayXXf::Zero(sequenceLength, states);
    m_decodingReturnValues.sumOverPairs11 = Eigen::ArrayXXf::Zero(sequenceLength, states);
  }
}

// Decodes all pairs. Returns a sum of all decoded posteriors (sequenceLength x states).
void HMM::decodeAll(int jobs, int jobInd)
{

  // auto t0 = std::chrono::high_resolution_clock().now();
  Timer timer;

  resetDecoding();

  const vector<Individual>& individuals = data.individuals;
  uint64 lastPercentage = -1;
  uint64 N = individuals.size();
  uint64 pairs = 0, pairsJob = 0;

  // calculate total number of pairs to decode
  uint64 totPairs;
  if (!decodingParams.withinOnly) {
    totPairs = 2 * N * N - N;
  } else {
    totPairs = N;
  }

  // figure out the range of pairs for this job number
  uint64 pairsStart = totPairs * (jobInd - 1) / jobs;
  uint64 pairsEnd = totPairs * jobInd / jobs;
  uint64 totPairsJob = pairsEnd - pairsStart;

  // alloc
  m_observationsBatch.clear();
  for (uint i = 0; i < individuals.size(); i++) {
    if (!decodingParams.withinOnly) {
      for (uint j = 0; j < i; j++) {
        // different individuals; decode 2 haps x 2 haps
        for (int iHap = 1; iHap <= 2; iHap++) {
          for (int jHap = 1; jHap <= 2; jHap++) {
            if (pairsStart <= pairs && pairs < pairsEnd) {
              PairObservations observations = makePairObs(individuals[i], iHap, individuals[j], jHap);
              if (noBatches) {
                decode(observations);
              } else {
                addToBatch(m_observationsBatch, observations);
              }
              pairsJob++;
            }
            pairs++;
            currPair = pairs;
          }
        }
      }
    }
    // this is the same individual; only decode across chromosomes
    if (pairsStart <= pairs && pairs < pairsEnd) {
      PairObservations observations = makePairObs(individuals[i], 1, individuals[i], 2);
      if (noBatches) {
        decode(observations);
      } else {
        addToBatch(m_observationsBatch, observations);
      }
      pairsJob++;
    }
    pairs++;
    currPair = pairs;
    uint64 percentage = 100 * pairsJob / totPairsJob;
    if (percentage != lastPercentage) {
      cout << "\rDecoding progress: " << percentage << "%"
           << "  (" << pairsJob << "/" << totPairsJob << ")" << flush;
    }
    lastPercentage = percentage;
  }

  // auto t1 = std::chrono::high_resolution_clock().now();
  // double ticksDecodeAll = t1 - t0;
  printf("\nDecoded %" PRIu64 " pairs in %.9f seconds.\n", pairsJob, timer.update_time());
  // print some stats (will remove)
  //    printPctTime("forward", ticksForward / ticksDecodeAll);
  //    printPctTime("backward", ticksBackward / ticksDecodeAll);
  //    printPctTime("combine", ticksCombine / ticksDecodeAll);
  //    printPctTime("sumOverPairs", ticksSumOverPairs / ticksDecodeAll);
  //    printPctTime("outputPerPair", ticksOutputPerPair / ticksDecodeAll);
  //    printPctTime("other", 1 - (ticksForward + ticksBackward + ticksCombine +
  //    ticksSumOverPairs + ticksOutputPerPair) / ticksDecodeAll);
  cout << flush;

  finishDecoding();
}

void HMM::decodePairs(const vector<uint>& individualsA, const vector<uint>& individualsB)
{
  if (individualsA.size() != individualsB.size()) {
    throw runtime_error("vector of A indicies must be the same size as vector of B indicies");
  }
  for (size_t i = 0; i < individualsA.size(); ++i) {
    decodePair(individualsA[i], individualsB[i]);
  }
}

void HMM::decodePair(const uint i, const uint j)
{
  const vector<Individual>& individuals = data.individuals;
  assert(i < individuals.size());
  assert(j < individuals.size());

  if (i != j) {
    // different individuals; decode 2 haps x 2 haps
    for (int iHap = 1; iHap <= 2; iHap++) {
      for (int jHap = 1; jHap <= 2; jHap++) {
        PairObservations observations = makePairObs(individuals[i], iHap, individuals[j], jHap);
        if (noBatches) {
          decode(observations);
        } else {
          addToBatch(m_observationsBatch, observations);
        }
      }
    }
  } else {
    // this is the same individual; only decode across chromosomes
    PairObservations observations = makePairObs(individuals[i], 1, individuals[i], 2);
    if (noBatches) {
      decode(observations);
    } else {
      addToBatch(m_observationsBatch, observations);
    }
  }
}

unsigned HMM::getFromPosition(unsigned from, const double cmDist)
{
  assert(cmDist > 0.0);

  double cumGenDist = 0.0;
  while (cumGenDist < cmDist && from > 0) {
    from--;
    cumGenDist += (data.geneticPositions[from + 1] - data.geneticPositions[from]) * 100.0;
  }
  return from;
}

unsigned HMM::getToPosition(unsigned to, const double cmDist)
{
  assert(cmDist > 0.0);

  double cumGenDist = 0.0;
  while (cumGenDist < cmDist && to < sequenceLength) {
    cumGenDist += (data.geneticPositions[to] - data.geneticPositions[to - 1]) * 100.0;
    to++;
  }
  return to;
}



uint HMM::getStateThreshold()
{
  uint result = 0u;
  const vector <float>& disc = m_decodingQuant.discretization;

  while (disc[result] < static_cast<float>(decodingParams.time) && result < m_decodingQuant.states) {
    result++;
  }
  return stateThreshold;
}

void HMM::finishDecoding()
{
  runLastBatch(m_observationsBatch);
  if (decodingParams.doPerPairPosteriorMean) {
    foutPosteriorMeanPerPair.close();
  }
  if (decodingParams.doPerPairMAP) {
    foutMAPPerPair.close();
  }
}

// add pair to batch and run if we have enough
void HMM::addToBatch(vector<PairObservations>& obsBatch, const PairObservations& observations)
{
  obsBatch.push_back(observations);
  if ((int)obsBatch.size() == m_batchSize) {
    // decodeBatch saves posteriors into m_alphaBuffer [sequenceLength x states x
    // m_batchSize]
    decodeBatch(obsBatch);
    augmentSumOverPairs(obsBatch, m_batchSize, m_batchSize);
    if (decodingParams.doPerPairMAP || decodingParams.doPerPairPosteriorMean) {
      writePerPairOutput(m_batchSize, m_batchSize, obsBatch);
    }
    obsBatch.clear();
  }
}

//// complete with leftover pairs
//void runLastBatchFastSMC(vector<PairObservations>& obsBatch)
//{
//  if (!obsBatch.empty()) {
//
//    int actualBatchSize = obsBatch.size();
//
//    // taking the maximum To position and the minimum From position in the batch
//    startBatch = *std::min_element(fromBatch.begin(), fromBatch.begin() + actualBatchSize);
//    endBatch = *std::max_element(toBatch.begin(), toBatch.begin() + actualBatchSize);
//
//    unsigned int from = getFromPosition(startBatch);
//    unsigned int to = getToPosition(endBatch);
//
//    for (uint i = 0; i < obsBatch.size(); i++) {
//      makeBits(obsBatch[i], from, to);
//    }
//
//    while (obsBatch.size() % VECX != 0) { // fill to size divisible by VECX
//      obsBatch.push_back(obsBatch.back());
//    }
//
//    int paddedBatchSize = obsBatch.size();
//
//    // decodeBatch saves posteriors into alphaBuffer [sequenceLength x states x paddedBatchSize]
//    decodeBatchFastSMC(obsBatch, actualBatchSize, paddedBatchSize, from, to);
//    obsBatch.clear();
//  }
//}

// complete with leftover pairs
void HMM::runLastBatch(vector<PairObservations>& obsBatch)
{
  if (obsBatch.empty())
    return;
  int actualBatchSize = static_cast<int>(obsBatch.size());
  while (obsBatch.size() % VECX != 0) // fill to size divisible by VECX
    obsBatch.push_back(obsBatch.back());
  int paddedBatchSize = static_cast<int>(obsBatch.size());
  // decodeBatch saves posteriors into m_alphaBuffer [sequenceLength x states x
  // paddedBatchSize]
  decodeBatch(obsBatch);
  augmentSumOverPairs(obsBatch, actualBatchSize, paddedBatchSize);
  if (decodingParams.doPerPairMAP || decodingParams.doPerPairPosteriorMean) {
    writePerPairOutput(actualBatchSize, paddedBatchSize, obsBatch);
  }
  obsBatch.clear();
}

// decode a batch
void HMM::decodeBatch(const vector<PairObservations>& obsBatch)
{

  int curBatchSize = static_cast<int>(obsBatch.size());

  float* obsIsZeroBatch = ALIGNED_MALLOC_FLOATS(sequenceLength * curBatchSize);
  float* obsIsTwoBatch = ALIGNED_MALLOC_FLOATS(sequenceLength * curBatchSize);

  for (long int pos = 0; pos < sequenceLength; pos++) {
    for (int v = 0; v < curBatchSize; v++) {
      obsIsZeroBatch[pos * curBatchSize + v] = (!obsBatch[v].obsBits[pos] ? 1.0f : 0.0f);
      obsIsTwoBatch[pos * curBatchSize + v] = (obsBatch[v].homMinorBits[pos] ? 1.0f : 0.0f);
    }
  }

  auto t0 = std::chrono::high_resolution_clock().now();

  // run forward
  forwardBatch(obsIsZeroBatch, obsIsTwoBatch, curBatchSize);

  auto t1 = std::chrono::high_resolution_clock().now();
  ticksForward += t1 - t0;

  // run backward
  backwardBatch(obsIsZeroBatch, obsIsTwoBatch, curBatchSize);

  auto t2 = std::chrono::high_resolution_clock().now();
  ticksBackward += t2 - t1;

  // combine (alpha * beta), normalize and store
  float* scale = obsIsZeroBatch; // reuse buffer but rename to be less confusing
  memset(scale, 0, sequenceLength * curBatchSize * sizeof(scale[0]));
#ifdef NO_SSE
  for (long int pos = 0; pos < sequenceLength; pos++) {
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v++) {
        long int ind = (pos * states + k) * curBatchSize + v;
        m_alphaBuffer[ind] *= m_betaBuffer[ind];
        scale[pos * curBatchSize + v] += m_alphaBuffer[ind];
      }
    }
  }
  for (long int pos = 0; pos < sequenceLength; pos++) {
    for (int v = 0; v < curBatchSize; v++) {
      scale[pos * curBatchSize + v] = 1.0f / scale[pos * curBatchSize + v];
    }
  }
  for (long int pos = 0; pos < sequenceLength; pos++) {
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v++) {
        m_alphaBuffer[(pos * states + k) * curBatchSize + v] *= scale[pos * curBatchSize + v];
      }
    }
  }
#else
  for (long int pos = 0; pos < sequenceLength; pos++) {
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v += VECX) {
        long int ind = (pos * states + k) * curBatchSize + v;
        FLOAT prod = MULT(LOAD(&m_alphaBuffer[ind]), LOAD(&m_betaBuffer[ind]));
        STORE(&m_alphaBuffer[ind], prod);
        STORE(&scale[pos * curBatchSize + v], ADD(LOAD(&scale[pos * curBatchSize + v]), prod));
      }
    }
  }
  for (long int pos = 0; pos < sequenceLength; pos++) {
    for (int v = 0; v < curBatchSize; v += VECX) {
      long int ind = pos * curBatchSize + v;
      STORE(&scale[ind], RECIPROCAL(LOAD(&scale[ind])));
    }
  }
  for (long int pos = 0; pos < sequenceLength; pos++) {
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v += VECX) {
        long int ind = (pos * states + k) * curBatchSize + v;
        STORE(&m_alphaBuffer[ind], MULT(LOAD(&m_alphaBuffer[ind]), LOAD(&scale[pos * curBatchSize + v])));
      }
    }
  }
#endif

  auto t3 = std::chrono::high_resolution_clock().now();
  ticksCombine += t3 - t2;
  ALIGNED_FREE(obsIsZeroBatch);
  ALIGNED_FREE(obsIsTwoBatch);
}

// compute scaling factor for an alpha vector
void HMM::scaleBatch(float* alpha, float* scalings, float* sums, int curBatchSize)
{
#ifdef NO_SSE
  // compute scaling (sum of current alpha vector)
  for (int k = 0; k < states; k++) {
    for (int v = 0; v < curBatchSize; v++) {
      sums[v] += alpha[k * curBatchSize + v];
    }
  }
  for (int v = 0; v < curBatchSize; v++) {
    scalings[v] = 1.0f / sums[v];
  }
#else
  // compute scaling (sum of current alpha vector)
  for (int k = 0; k < states; k++)
    for (int v = 0; v < curBatchSize; v += VECX) {
      STORE(&sums[v], ADD(LOAD(&sums[v]), LOAD(&alpha[k * curBatchSize + v])));
    }
  for (int v = 0; v < curBatchSize; v++)
    scalings[v] = 1.0f / sums[v];
#endif
}

// apply scaling factor to alpha/beta vector
void HMM::applyScaling(float* vec, float* scalings, int curBatchSize)
{
#ifdef NO_SSE
  // normalize current alpha vector to 1
  for (int k = 0; k < states; k++) {
    for (int v = 0; v < curBatchSize; v++) {
      vec[k * curBatchSize + v] *= scalings[v];
    }
  }
#else
  // normalize current alpha vector to 1
  for (int k = 0; k < states; k++)
    for (int v = 0; v < curBatchSize; v += VECX) {
      STORE(&vec[k * curBatchSize + v], MULT(LOAD(&vec[k * curBatchSize + v]), LOAD(&scalings[v])));
    }
#endif
}

// forward step
void HMM::forwardBatch(const float* obsIsZeroBatch, const float* obsIsTwoBatch, int curBatchSize)
{

  assert(curBatchSize % VECX == 0);

  float* alphaC = ALIGNED_MALLOC_FLOATS(states * curBatchSize);
  float* AU = ALIGNED_MALLOC_FLOATS(curBatchSize);

  // fill pos=0 in alpha
  for (int k = 0; k < states; k++) {
    for (int v = 0; v < curBatchSize; v++) {
      float firstEmission = emission1AtSite[0][k] + emission0minus1AtSite[0][k] * obsIsZeroBatch[v] +
                            emission2minus0AtSite[0][k] * obsIsTwoBatch[v];
      m_alphaBuffer[k * curBatchSize + v] = m_decodingQuant.initialStateProb[k] * firstEmission;
    }
  }

  float* sums = AU; // reuse buffer but rename to be less confusing
  memset(sums, 0, curBatchSize * sizeof(sums[0]));
  scaleBatch(m_alphaBuffer, m_scalingBuffer, sums, curBatchSize);
  applyScaling(m_alphaBuffer, m_scalingBuffer, curBatchSize);

  // Induction Step:
  float lastGeneticPos = data.geneticPositions[0];
  int lastPhysicalPos = data.physicalPositions[0];

  for (long int pos = 1; pos < sequenceLength; pos++) {
    // get distances and rates
    float recDistFromPrevious = roundMorgans(std::max(minGenetic, data.geneticPositions[pos] - lastGeneticPos));
    float currentRecRate = roundMorgans(data.recRateAtMarker[pos]);
    float* previousAlpha = &m_alphaBuffer[(pos - 1) * states * curBatchSize];
    float* nextAlpha = &m_alphaBuffer[pos * states * curBatchSize];
    if (decodingParams.decodingSequence) {
      int physDistFromPreviousMinusOne = roundPhysical(data.physicalPositions[pos] - lastPhysicalPos - 1);
      float recDistFromPreviousMinusOne = roundMorgans(std::max(minGenetic, recDistFromPrevious - currentRecRate));
      vector<float> homozEmission = m_decodingQuant.homozygousEmissionMap.at(physDistFromPreviousMinusOne);
      getNextAlphaBatched(recDistFromPreviousMinusOne, alphaC, curBatchSize, previousAlpha, pos, m_allZeros, m_allZeros,
                          AU, nextAlpha, homozEmission, homozEmission, homozEmission);
      previousAlpha = nextAlpha;
      getNextAlphaBatched(currentRecRate, alphaC, curBatchSize, previousAlpha, pos, obsIsZeroBatch, obsIsTwoBatch, AU,
                          nextAlpha, emission1AtSite[pos], emission0minus1AtSite[pos], emission2minus0AtSite[pos]);
    } else {
      getNextAlphaBatched(recDistFromPrevious, alphaC, curBatchSize, previousAlpha, pos, obsIsZeroBatch, obsIsTwoBatch,
                          AU, nextAlpha, emission1AtSite[pos], emission0minus1AtSite[pos], emission2minus0AtSite[pos]);
    }
    float* sums = AU; // reuse buffer but rename to be less confusing
    memset(sums, 0, curBatchSize * sizeof(sums[0]));
    if (pos % scalingSkip == 0) {
      scaleBatch(nextAlpha, m_scalingBuffer, sums, curBatchSize);
      applyScaling(nextAlpha, m_scalingBuffer, curBatchSize);
    }
    // update distances
    lastGeneticPos = data.geneticPositions[pos];
    lastPhysicalPos = data.physicalPositions[pos];
  }

  ALIGNED_FREE(AU);
  ALIGNED_FREE(alphaC);
}

// compute next alpha vector in linear time
void HMM::getNextAlphaBatched(float recDistFromPrevious, float* alphaC, int curBatchSize, const float* previousAlpha,
                              uint pos, const float* obsIsZeroBatch, const float* obsIsTwoBatch, float* AU,
                              float* nextAlpha, const vector<float>& emission1AtSite,
                              const vector<float>& emission0minus1AtSite, const vector<float>& emission2minus0AtSite)
{

  const float* B = &m_decodingQuant.Bvectors.at(recDistFromPrevious)[0];
  const float* U = &m_decodingQuant.Uvectors.at(recDistFromPrevious)[0];
  const float* D = &m_decodingQuant.Dvectors.at(recDistFromPrevious)[0];

  memcpy(&alphaC[(states - 1) * curBatchSize], &previousAlpha[(states - 1) * curBatchSize],
         curBatchSize * sizeof(alphaC[0]));

  for (int k = states - 2; k >= 0; k--) {
#ifdef NO_SSE
    for (int v = 0; v < curBatchSize; v++) {
      alphaC[k * curBatchSize + v] = alphaC[(k + 1) * curBatchSize + v] + previousAlpha[k * curBatchSize + v];
    }
#else
    for (int v = 0; v < curBatchSize; v += VECX) {
      FLOAT alphaC_kp1_v = LOAD(&alphaC[(k + 1) * curBatchSize + v]);
      FLOAT prevAlpha_k_v = LOAD(&previousAlpha[k * curBatchSize + v]);
      STORE(&alphaC[k * curBatchSize + v], ADD(alphaC_kp1_v, prevAlpha_k_v));
    }
#endif
  }

  memset(AU, 0, curBatchSize * sizeof(AU[0]));
  for (int k = 0; k < states; k++) {

#ifdef NO_SSE
    for (int v = 0; v < curBatchSize; v++) {
      if (k)
        AU[v] = U[k - 1] * previousAlpha[(k - 1) * curBatchSize + v] + m_decodingQuant.columnRatios[k - 1] * AU[v];
      float term = AU[v] + D[k] * previousAlpha[k * curBatchSize + v];
      if (k < states - 1) {
        term += B[k] * alphaC[(k + 1) * curBatchSize + v];
      }
      float currentEmission_k = emission1AtSite[k] + emission0minus1AtSite[k] * obsIsZeroBatch[pos * curBatchSize + v] +
                                emission2minus0AtSite[k] * obsIsTwoBatch[pos * curBatchSize + v];
      nextAlpha[k * curBatchSize + v] = currentEmission_k * term;
    }
#else
#ifdef AVX512
    FLOAT D_k = LOAD1(D[k]);
    FLOAT B_k = LOAD1(B[k]);
    FLOAT em1Prob_k = LOAD1(emission1AtSite[k]);
    FLOAT em0minus1Prob_k = LOAD1(emission0minus1AtSite[k]);
    FLOAT em2minus0Prob_k = LOAD1(emission2minus0AtSite[k]);
#else
    FLOAT D_k = LOAD1(&D[k]);
    FLOAT B_k = LOAD1(&B[k]);
    FLOAT em1Prob_k = LOAD1(&emission1AtSite[k]);
    FLOAT em0minus1Prob_k = LOAD1(&emission0minus1AtSite[k]);
    FLOAT em2minus0Prob_k = LOAD1(&emission2minus0AtSite[k]);
#endif

    FLOAT Ukm1, colRatios_km1;
    if (k) {
#ifdef AVX512
      Ukm1 = LOAD1(U[k - 1]);
      colRatios_km1 = LOAD1(m_decodingQuant.columnRatios[k - 1]);
#else
      Ukm1 = LOAD1(&U[k - 1]);
      colRatios_km1 = LOAD1(&m_decodingQuant.columnRatios[k - 1]);
#endif
    }
    for (int v = 0; v < curBatchSize; v += VECX) {
      FLOAT AU_v;
      if (k) {
        FLOAT term1 = MULT(Ukm1, LOAD(&previousAlpha[(k - 1) * curBatchSize + v]));
        FLOAT term2 = MULT(colRatios_km1, LOAD(&AU[v]));
        AU_v = ADD(term1, term2);
        STORE(&AU[v], AU_v);
      } else
        AU_v = LOAD(&AU[v]);

      FLOAT term = ADD(AU_v, MULT(D_k, LOAD(&previousAlpha[k * curBatchSize + v])));
      if (k < states - 1) { // TODO: just extend B and alphaC?
        term = ADD(term, MULT(B_k, LOAD(&alphaC[(k + 1) * curBatchSize + v])));
      }
      FLOAT currentEmission_k =
          ADD(ADD(em1Prob_k, MULT(em0minus1Prob_k, LOAD(&obsIsZeroBatch[pos * curBatchSize + v]))),
              MULT(em2minus0Prob_k, LOAD(&obsIsTwoBatch[pos * curBatchSize + v])));
      if (v == 0) {
      }
      STORE(&nextAlpha[k * curBatchSize + v], MULT(currentEmission_k, term));
    }
#endif
  }
}

// backward step
void HMM::backwardBatch(const float* obsIsZeroBatch, const float* obsIsTwoBatch, int curBatchSize)
{

  // fill pos=sequenceLenght-1 in beta
  for (int k = 0; k < states; k++) {
    for (int v = 0; v < curBatchSize; v++) {
      m_betaBuffer[((sequenceLength - 1) * states + k) * curBatchSize + v] = 1.0f;
    }
  }
  float* sums = ALIGNED_MALLOC_FLOATS(curBatchSize);
  memset(sums, 0, curBatchSize * sizeof(sums[0]));
  scaleBatch(m_betaBuffer, m_scalingBuffer, sums, curBatchSize);
  applyScaling(m_betaBuffer, m_scalingBuffer, curBatchSize);

  // Induction Step:
  float* BL = ALIGNED_MALLOC_FLOATS(curBatchSize);
  float* BU = ALIGNED_MALLOC_FLOATS(states * curBatchSize);
  float* vec = ALIGNED_MALLOC_FLOATS(states * curBatchSize);

  float lastGeneticPos = data.geneticPositions[sequenceLength - 1];
  int lastPhysicalPos = data.physicalPositions[sequenceLength - 1];

  for (long int pos = sequenceLength - 2; pos >= 0; pos--) {
    // get distances and rates
    float recDistFromPrevious = roundMorgans(std::max(minGenetic, lastGeneticPos - data.geneticPositions[pos]));
    float currentRecRate = roundMorgans(data.recRateAtMarker[pos]);
    float* currentBeta = &m_betaBuffer[pos * states * curBatchSize];
    float* lastComputedBeta = &m_betaBuffer[(pos + 1) * states * curBatchSize];
    if (decodingParams.decodingSequence) {
      int physDistFromPreviousMinusOne = roundPhysical(lastPhysicalPos - data.physicalPositions[pos] - 1);
      float recDistFromPreviousMinusOne = roundMorgans(std::max(minGenetic, recDistFromPrevious - currentRecRate));
      vector<float> homozEmission = m_decodingQuant.homozygousEmissionMap.at(physDistFromPreviousMinusOne);
      getPreviousBetaBatched(recDistFromPreviousMinusOne, curBatchSize, lastComputedBeta, pos, m_allZeros, m_allZeros,
                             vec, BU, BL, currentBeta, homozEmission, homozEmission, homozEmission);
      lastComputedBeta = currentBeta;
      getPreviousBetaBatched(currentRecRate, curBatchSize, lastComputedBeta, pos, obsIsZeroBatch, obsIsTwoBatch, vec,
                             BU, BL, currentBeta, emission1AtSite[pos + 1], emission0minus1AtSite[pos + 1],
                             emission2minus0AtSite[pos + 1]);
    } else {
      getPreviousBetaBatched(recDistFromPrevious, curBatchSize, lastComputedBeta, pos, obsIsZeroBatch, obsIsTwoBatch,
                             vec, BU, BL, currentBeta, emission1AtSite[pos + 1], emission0minus1AtSite[pos + 1],
                             emission2minus0AtSite[pos + 1]);
    }
    if (pos % scalingSkip == 0) {
      // normalize betas using alpha scaling
      memset(sums, 0, curBatchSize * sizeof(sums[0]));
      scaleBatch(currentBeta, m_scalingBuffer, sums, curBatchSize);
      applyScaling(currentBeta, m_scalingBuffer, curBatchSize);
    }
    // update distances
    lastGeneticPos = data.geneticPositions[pos];
    lastPhysicalPos = data.physicalPositions[pos];
  }

  ALIGNED_FREE(vec);
  ALIGNED_FREE(BU);
  ALIGNED_FREE(BL);
}

// compute previous beta vector in linear time
void HMM::getPreviousBetaBatched(float recDistFromPrevious, int curBatchSize, const float* lastComputedBeta, int pos,
                                 const float* obsIsZeroBatch, const float* obsIsTwoBatch, float* vec, float* BU,
                                 float* BL, float* currentBeta, const vector<float>& emission1AtSite,
                                 const vector<float>& emission0minus1AtSite, const vector<float>& emission2minus0AtSite)
{
  const vector<float>& B = m_decodingQuant.Bvectors.at(recDistFromPrevious);
  const vector<float>& U = m_decodingQuant.Uvectors.at(recDistFromPrevious);
  const vector<float>& RR = m_decodingQuant.rowRatioVectors.at(recDistFromPrevious);
  const vector<float>& D = m_decodingQuant.Dvectors.at(recDistFromPrevious);
#ifdef NO_SSE

  for (int k = 0; k < states; k++) {
    for (int v = 0; v < curBatchSize; v++) {
      float currentEmission_k = emission1AtSite[k] +
                                emission0minus1AtSite[k] * obsIsZeroBatch[(pos + 1) * curBatchSize + v] +
                                emission2minus0AtSite[k] * obsIsTwoBatch[(pos + 1) * curBatchSize + v];
      vec[k * curBatchSize + v] = lastComputedBeta[k * curBatchSize + v] * currentEmission_k;
    }
  }

#else
  for (int k = 0; k < states; k++) {
#ifdef AVX512
    FLOAT em1Prob_k = LOAD1(emission1AtSite[k]);
    FLOAT em0minus1Prob_k = LOAD1(emission0minus1AtSite[k]);
    FLOAT em2minus0Prob_k = LOAD1(emission2minus0AtSite[k]);
#else
    FLOAT em1Prob_k = LOAD1(&emission1AtSite[k]);
    FLOAT em0minus1Prob_k = LOAD1(&emission0minus1AtSite[k]);
    FLOAT em2minus0Prob_k = LOAD1(&emission2minus0AtSite[k]);
#endif
    for (int v = 0; v < curBatchSize; v += VECX) {
      FLOAT currentEmission_k =
          ADD(ADD(em1Prob_k, MULT(em0minus1Prob_k, LOAD(&obsIsZeroBatch[(pos + 1) * curBatchSize + v]))),
              MULT(em2minus0Prob_k, LOAD(&obsIsTwoBatch[(pos + 1) * curBatchSize + v])));
      STORE(&vec[k * curBatchSize + v], MULT(LOAD(&lastComputedBeta[k * curBatchSize + v]), currentEmission_k));
    }
  }
#endif

  memset(BU + (states - 1) * curBatchSize, 0, curBatchSize * sizeof(BU[0]));
#ifdef NO_SSE
  for (int k = states - 2; k >= 0; k--)
    for (int v = 0; v < curBatchSize; v++)
      BU[k * curBatchSize + v] = U[k] * vec[(k + 1) * curBatchSize + v] + RR[k] * BU[(k + 1) * curBatchSize + v];
#else
  for (int k = states - 2; k >= 0; k--) {
#ifdef AVX512
    FLOAT U_k = LOAD1(U[k]);
    FLOAT RR_k = LOAD1(RR[k]);
#else
    FLOAT U_k = LOAD1(&U[k]);
    FLOAT RR_k = LOAD1(&RR[k]);
#endif
    for (int v = 0; v < curBatchSize; v += VECX) {
      FLOAT term1 = MULT(U_k, LOAD(&vec[(k + 1) * curBatchSize + v]));
      FLOAT term2 = MULT(RR_k, LOAD(&BU[(k + 1) * curBatchSize + v]));
      STORE(&BU[k * curBatchSize + v], ADD(term1, term2));
    }
  }
#endif

  memset(BL, 0, curBatchSize * sizeof(BL[0]));
#ifdef NO_SSE
  for (int k = 0; k < states; k++) {
    for (int v = 0; v < curBatchSize; v++) {
      if (k)
        BL[v] += B[k - 1] * vec[(k - 1) * curBatchSize + v];
      currentBeta[k * curBatchSize + v] = BL[v] + D[k] * vec[k * curBatchSize + v] + BU[k * curBatchSize + v];
    }
  }
#else
  for (int k = 0; k < states; k++) {
#ifdef AVX512
    FLOAT D_k = LOAD1(D[k]);
    FLOAT B_km1;
    if (k)
      B_km1 = LOAD1(B[k - 1]);
#else
    FLOAT D_k = LOAD1(&D[k]);
    FLOAT B_km1;
    if (k)
      B_km1 = LOAD1(&B[k - 1]);
#endif
    for (int v = 0; v < curBatchSize; v += VECX) {
      FLOAT BL_v = LOAD(&BL[v]);
      if (k) {
        BL_v = ADD(BL_v, MULT(B_km1, LOAD(&vec[(k - 1) * curBatchSize + v])));
        STORE(&BL[v], BL_v);
      }
      STORE(&currentBeta[k * curBatchSize + v],
            ADD(BL_v, ADD(MULT(D_k, LOAD(&vec[k * curBatchSize + v])), LOAD(&BU[k * curBatchSize + v]))));
    }
  }
#endif
}

// --posteriorSums
void HMM::augmentSumOverPairs(vector<PairObservations>& obsBatch, int actualBatchSize, int paddedBatchSize)
{

  auto t0 = std::chrono::high_resolution_clock().now();

  if (!decodingParams.doPosteriorSums && !decodingParams.doMajorMinorPosteriorSums)
    return;

  for (long int pos = 0; pos < sequenceLength; pos++) {
    for (int k = 0; k < states; k++) {
      float sum = 0;
      float sum00 = 0;
      float sum01 = 0;
      float sum11 = 0;
      for (int v = 0; v < actualBatchSize; v++) { // only loop over actual (not padding) pairs
        float posterior_pos_state_pair = m_alphaBuffer[(pos * states + k) * paddedBatchSize + v];
        if (decodingParams.doPosteriorSums) {
          sum += posterior_pos_state_pair;
        }
        if (decodingParams.doMajorMinorPosteriorSums) {
          if (obsBatch[v].homMinorBits[pos] == 1)
            sum11 += posterior_pos_state_pair;
          else if (obsBatch[v].obsBits[pos] == 0)
            sum00 += posterior_pos_state_pair;
          else
            sum01 += posterior_pos_state_pair;
        }
      }
      if (decodingParams.doPosteriorSums) {
        m_decodingReturnValues.sumOverPairs(pos, k) += sum;
      }
      if (decodingParams.doMajorMinorPosteriorSums) {
        m_decodingReturnValues.sumOverPairs00(pos, k) += sum00;
        m_decodingReturnValues.sumOverPairs01(pos, k) += sum01;
        m_decodingReturnValues.sumOverPairs11(pos, k) += sum11;
      }
    }
  }

  auto t1 = std::chrono::high_resolution_clock().now();
  ticksSumOverPairs += t1 - t0;
}

// will eventually write binary output instead of gzipped
void HMM::writePerPairOutput(int actualBatchSize, int paddedBatchSize, const vector<PairObservations>& obsBatch)
{

  auto t0 = std::chrono::high_resolution_clock().now();

  if (decodingParams.doPerPairMAP) {
    memset(MAP, 0, sequenceLength * actualBatchSize * sizeof(MAP[0]));
  }
  if (decodingParams.doPerPairPosteriorMean) {
    memset(meanPost, 0, sequenceLength * actualBatchSize * sizeof(meanPost[0]));
  }

  for (long int pos = 0; pos < sequenceLength; pos++) {
    if (decodingParams.doPerPairMAP) {
      memset(currentMAPValue, 0, actualBatchSize * sizeof(currentMAPValue[0]));
    }
    for (ushort k = 0; k < states; k++) {
      for (int v = 0; v < actualBatchSize; v++) {
        float posterior_pos_state_pair = m_alphaBuffer[(pos * states + k) * paddedBatchSize + v];
        if (decodingParams.doPerPairPosteriorMean) {
          meanPost[pos * actualBatchSize + v] += posterior_pos_state_pair * expectedCoalTimes[k];
        }
        if (decodingParams.doPerPairMAP && currentMAPValue[v] < posterior_pos_state_pair) {
          MAP[pos * actualBatchSize + v] = k;
          currentMAPValue[v] = posterior_pos_state_pair;
        }
      }
    }
  }

  if (decodingParams.doPerPairPosteriorMean) {
    for (int v = 0; v < actualBatchSize; v++) {
//      foutPosteriorMeanPerPair << obsBatch[v].iName << " " << obsBatch[v].iHap << " ";
//      foutPosteriorMeanPerPair << obsBatch[v].jName << " " << obsBatch[v].jHap;
      for (long int pos = 0; pos < sequenceLength; pos++) {
        foutPosteriorMeanPerPair << " " << meanPost[pos * actualBatchSize + v];
      }
      foutPosteriorMeanPerPair << endl;
    }
  }

  if (decodingParams.doPerPairMAP) {
    for (int v = 0; v < actualBatchSize; v++) {
//      foutMAPPerPair << obsBatch[v].iName << " " << obsBatch[v].iHap << " ";
//      foutMAPPerPair << obsBatch[v].jName << " " << obsBatch[v].jHap;
      for (long int pos = 0; pos < sequenceLength; pos++) {
        foutMAPPerPair << " " << MAP[pos * actualBatchSize + v];
      }
      foutMAPPerPair << endl;
    }
  }

  auto t1 = std::chrono::high_resolution_clock().now();
  ticksOutputPerPair += t1 - t0;
}

// *********************************************************************
// non-batched computations (for debugging and pedagogical reasons only)
// *********************************************************************

float HMM::printVector(const vector<float>& vec)
{
  float sum = 0.f;
  for (uint i = 0; i < vec.size(); i++) {
    cout << vec[i] << "\t";
  }
  cout << endl;
  return sum;
}

float HMM::getSumOfVector(const vector<float>& vec)
{
  float sum = 0.f;
  for (uint i = 0; i < vec.size(); i++) {
    sum += vec[i];
  }
  return sum;
}

vector<float> HMM::elementWiseMultVectorScalar(const vector<float>& vec, float val)
{
  vector<float> ret(vec.size());
  for (uint i = 0; i < vec.size(); i++) {
    ret[i] = vec[i] * val;
  }
  return ret;
}

vector<float> HMM::elementWiseMultVectorVector(const vector<float>& vec, const vector<float>& factors)
{
  vector<float> ret(vec.size());
  for (uint i = 0; i < vec.size(); i++) {
    ret[i] = vec[i] * factors[i];
  }
  return ret;
}

vector<vector<float>> HMM::elementWiseMultMatrixMatrix(const vector<vector<float>>& matrix1,
                                                       const vector<vector<float>>& matrix2)
{
  vector<vector<float>> ret(matrix1.size(), vector<float>(matrix1[0].size()));
  for (uint i = 0; i < matrix1.size(); i++) {
    for (uint j = 0; j < matrix1[0].size(); j++) {
      ret[i][j] = matrix1[i][j] * matrix2[i][j];
    }
  }
  return ret;
}

vector<vector<float>> HMM::normalizeMatrixColumns(const vector<vector<float>>& matrix)
{
  vector<vector<float>> ret(matrix.size(), vector<float>(matrix[0].size()));
  for (uint j = 0; j < matrix[0].size(); j++) {
    float sum = 0.f;
    for (uint i = 0; i < matrix.size(); i++) {
      sum += matrix[i][j];
    }
    for (uint i = 0; i < matrix.size(); i++) {
      ret[i][j] = matrix[i][j] / sum;
    }
  }
  return ret;
}

void HMM::fillMatrixColumn(vector<vector<float>>& matrix, const vector<float>& vec, long int pos)
{
  for (uint i = 0; i < vec.size(); i++) {
    matrix[i][pos] = vec[i];
  }
}

vector<vector<float>> HMM::decode(const PairObservations& observations)
{

  auto t0 = std::chrono::high_resolution_clock().now();
  vector<vector<float>> forwardOut = forward(observations);
  auto t1 = std::chrono::high_resolution_clock().now();
  ticksForward += t1 - t0;

  vector<vector<float>> backwardOut = backward(observations);
  auto t2 = std::chrono::high_resolution_clock().now();
  ticksBackward += t2 - t1;

  vector<vector<float>> posterior = elementWiseMultMatrixMatrix(forwardOut, backwardOut);
  posterior = normalizeMatrixColumns(posterior);

  auto t3 = std::chrono::high_resolution_clock().now();
  ticksCombine += t3 - t2;

  if (decodingParams.doPosteriorSums) {
    for (uint k = 0; k < m_decodingQuant.states; k++) {
      for (long int pos = 0; pos < sequenceLength; pos++) {
        m_decodingReturnValues.sumOverPairs(pos, k) += posterior[k][pos];
      }
    }
  }
  return posterior;
}

// Instead of returning a whole posterior, return the MAP and posterior mean
pair<vector<float>, vector<float>> HMM::decodeSummarize(const PairObservations& observations)
{
  vector<vector<float>> posterior = decode(observations);
  size_t num_discretizations = m_decodingQuant.expectedTimes.size();
  assert(num_discretizations == posterior.size());

  vector<float> posterior_map(posterior[0].size(), 0);
  vector<float> posterior_max_so_far(posterior[0].size(), 0);
  vector<float> posterior_mean(posterior[0].size(), 0);
  for (int i = 0; i < posterior.size(); ++i) {
    for (int j = 0; j < posterior[0].size(); ++j) {
      posterior_mean[j] += posterior[i][j] * m_decodingQuant.expectedTimes[i];
      if (posterior[i][j] > posterior_max_so_far[j]) {
        posterior_max_so_far[j] = posterior[i][j];
        posterior_map[j] = m_decodingQuant.expectedTimes[i];
      }
    }
  }
  return std::make_pair(posterior_map, posterior_mean);
}

float HMM::roundMorgans(float gen)
{
  float gene1e10 = gen * 1e10f;
  int L10 = std::max(0, (int)floor(log10(gene1e10)) - precision);
  float factor = static_cast<float>(pow(10, L10));
  float rounded = round(gene1e10 / factor) * factor;
  return std::max(minGenetic, rounded / 1e10f);
}

int HMM::roundPhysical(int phys)
{
  // Since HMM for sequence uses distance-1, it can be -1.
  if (phys < -1) {
    cerr << "ERROR. Int overflow " << phys << endl;
    exit(1);
  }
  // map 0 or -1 to 1. Since HMM for sequence uses distance-1, it can be -1.
  phys = std::max(1, phys);
  int L10 = std::max(0, (int)floor(log10(phys)) - precision);
  int factor = (int)pow(10, L10);
  int rounded = static_cast<int>(round(phys / (float)factor)) * factor;
  return rounded;
}

vector<float> HMM::getEmission(int pos, int distinguished, int undistinguished, int emissionIndex)
{
  vector<float> emission;
  if (!useCSFSatThisPosition[pos]) {
    // this position is not a CSFS position, use compressed or classic
    emission = decodingParams.decodingSequence ? m_decodingQuant.classicEmissionTable[emissionIndex]
                                               : m_decodingQuant.compressedEmissionTable[emissionIndex];
  } else {
    // this position is a CSFS position
    if (decodingParams.foldData) {
      emission = decodingParams.decodingSequence
                     ? m_decodingQuant.foldedCSFSmap[undistinguished][emissionIndex]
                     : m_decodingQuant.foldedAscertainedCSFSmap[undistinguished][emissionIndex];
    } else {
      emission = decodingParams.decodingSequence ? m_decodingQuant.CSFSmap[undistinguished][distinguished]
                                                 : m_decodingQuant.ascertainedCSFSmap[undistinguished][distinguished];
    }
  }
  return emission;
}

// forward step
vector<vector<float>> HMM::forward(const PairObservations& observations)
{

  vector<vector<float>> alpha(states, vector<float>(sequenceLength));

  uint emissionIndex = observations.obsBits[0] ? 1 : 0;
  // if both samples are carriers, there are two distinguished, otherwise, it's the or
  // (use previous xor). This affects the number of undistinguished for the site
  uint distinguished = observations.homMinorBits[0] ? 2 : emissionIndex;
  uint undistinguished = useCSFSatThisPosition[0] ? data.undistinguishedCounts[0][distinguished] : -1;
  vector<float> emission = getEmission(0, distinguished, undistinguished, emissionIndex);

  vector<float> firstAlpha = elementWiseMultVectorVector(m_decodingQuant.initialStateProb, emission);

  // cumpute scaling (sum of current alpha vector)
  float m_scalingBuffer = 1.f / getSumOfVector(firstAlpha);
  // normalize current alpha vector to 1
  firstAlpha = elementWiseMultVectorScalar(firstAlpha, m_scalingBuffer);

  fillMatrixColumn(alpha, firstAlpha, 0);
  // Induction Step:
  vector<float> nextAlpha(states);
  vector<float> alphaC(states + 1);
  vector<float> previousAlpha = firstAlpha;
  float lastGeneticPos = data.geneticPositions[0];
  int lastPhysicalPos = data.physicalPositions[0];

  for (long int pos = 1; pos < sequenceLength; pos++) {
    auto t0 = std::chrono::high_resolution_clock().now();
    // get distances and rates
    float recDistFromPrevious = roundMorgans(std::max(minGenetic, data.geneticPositions[pos] - lastGeneticPos));
    float currentRecRate = roundMorgans(data.recRateAtMarker[pos]);
    // if both samples are carriers, there are two distinguished, otherwise, it's the
    // or (use previous xor). This affects the number of undistinguished for the site
    float obsIsZero = !observations.obsBits[pos] ? 1.0f : 0.0f;
    float obsIsHomMinor = observations.homMinorBits[pos] ? 1.0f : 0.0f;

    if (decodingParams.decodingSequence) {
      int physDistFromPreviousMinusOne = roundPhysical(data.physicalPositions[pos] - lastPhysicalPos - 1);
      float recDistFromPreviousMinusOne = roundMorgans(std::max(minGenetic, recDistFromPrevious - currentRecRate));
      vector<float> homozEmission = m_decodingQuant.homozygousEmissionMap.at(physDistFromPreviousMinusOne);
      getNextAlpha(recDistFromPreviousMinusOne, alphaC, previousAlpha, nextAlpha, homozEmission, homozEmission,
                   homozEmission, 0.0f, 0.0f);
      previousAlpha = nextAlpha;
      getNextAlpha(currentRecRate, alphaC, previousAlpha, nextAlpha, emission1AtSite[pos], emission0minus1AtSite[pos],
                   emission2minus0AtSite[pos], obsIsZero, obsIsHomMinor);
    } else {
      getNextAlpha(recDistFromPrevious, alphaC, previousAlpha, nextAlpha, emission1AtSite[pos],
                   emission0minus1AtSite[pos], emission2minus0AtSite[pos], obsIsZero, obsIsHomMinor);
    }
    auto t1 = std::chrono::high_resolution_clock().now();
    if (pos % scalingSkip == 0) {
      // compute scaling (sum of current alpha vector)
      m_scalingBuffer = 1.f / getSumOfVector(nextAlpha);
      // normalize current alpha vector to 1
      nextAlpha = elementWiseMultVectorScalar(nextAlpha, m_scalingBuffer);
    }
    fillMatrixColumn(alpha, nextAlpha, pos);
    previousAlpha = nextAlpha;
    auto t2 = std::chrono::high_resolution_clock().now();
    t1sum += t1 - t0;
    t2sum += t2 - t1;
    // update distances
    lastGeneticPos = data.geneticPositions[pos];
    lastPhysicalPos = data.physicalPositions[pos];
  }
  return alpha;
}

void HMM::getNextAlpha(float recDistFromPrevious, vector<float>& alphaC, vector<float>& previousAlpha,
                       vector<float>& nextAlpha, vector<float>& emission1AtSite, vector<float>& emission0minus1AtSite,
                       vector<float>& emission2minus0AtSite, float obsIsZero, float obsIsHomMinor)
{
  alphaC[m_decodingQuant.states - 1] = previousAlpha[m_decodingQuant.states - 1];
  for (int k = m_decodingQuant.states - 2; k >= 0; k--) {
    alphaC[k] = alphaC[k + 1] + previousAlpha[k];
  }
  const float* B = &m_decodingQuant.Bvectors.at(recDistFromPrevious)[0];
  const float* U = &m_decodingQuant.Uvectors.at(recDistFromPrevious)[0];
  const float* D = &m_decodingQuant.Dvectors.at(recDistFromPrevious)[0];
  float AUc = 0;
  for (uint k = 0; k < m_decodingQuant.states; k++) {
    if (k)
      AUc = U[k - 1] * previousAlpha[k - 1] + m_decodingQuant.columnRatios[k - 1] * AUc;
    float term = AUc + D[k] * previousAlpha[k];
    if (k < m_decodingQuant.states - 1)
      term += B[k] * alphaC[k + 1];
    float currentEmission_k =
        emission1AtSite[k] + emission0minus1AtSite[k] * obsIsZero + emission2minus0AtSite[k] * obsIsHomMinor;
    nextAlpha[k] = currentEmission_k * term;
  }
}

// backward step
vector<vector<float>> HMM::backward(const PairObservations& observations)
{

  vector<vector<float>> beta(states, vector<float>(sequenceLength));

  vector<float> lastBeta(states);
  for (uint i = 0; i < lastBeta.size(); i++) {
    lastBeta[i] = 1.f;
  }
  // normalize current alpha vector to 1
  float m_scalingBuffer = 1.f / getSumOfVector(lastBeta);
  lastBeta = elementWiseMultVectorScalar(lastBeta, m_scalingBuffer);
  fillMatrixColumn(beta, lastBeta, sequenceLength - 1);
  // Induction Step:
  vector<float> currentBeta(states);
  vector<float> BL(states);
  vector<float> BU(states);
  vector<float> lastComputedBeta = lastBeta;
  float lastGeneticPos = data.geneticPositions[sequenceLength - 1];
  int lastPhysicalPos = data.physicalPositions[sequenceLength - 1];
  for (long int pos = sequenceLength - 2; pos >= 0; pos--) {
    // get distances and rates
    float recDistFromPrevious = roundMorgans(std::max(minGenetic, lastGeneticPos - data.geneticPositions[pos]));
    float currentRecRate = roundMorgans(data.recRateAtMarker[pos]);
    // if both samples are carriers, there are two distinguished, otherwise, it's the
    // or (use previous xor). This affects the number of undistinguished for the site
    float obsIsZero = !observations.obsBits[pos + 1] ? 1.0f : 0.0f;
    float obsIsHomMinor = observations.homMinorBits[pos + 1] ? 1.0f : 0.0f;
    if (decodingParams.decodingSequence) {
      int physDistFromPreviousMinusOne = roundPhysical(lastPhysicalPos - data.physicalPositions[pos] - 1);
      float recDistFromPreviousMinusOne = roundMorgans(std::max(minGenetic, recDistFromPrevious - currentRecRate));
      vector<float> homozEmission = m_decodingQuant.homozygousEmissionMap.at(physDistFromPreviousMinusOne);
      getPreviousBeta(recDistFromPreviousMinusOne, lastComputedBeta, BL, BU, currentBeta, homozEmission, homozEmission,
                      homozEmission, 0.0f, 0.0f);
      lastComputedBeta = currentBeta;
      getPreviousBeta(currentRecRate, lastComputedBeta, BL, BU, currentBeta, emission1AtSite[pos + 1],
                      emission0minus1AtSite[pos + 1], emission2minus0AtSite[pos + 1], obsIsZero, obsIsHomMinor);
    } else {
      getPreviousBeta(recDistFromPrevious, lastComputedBeta, BL, BU, currentBeta, emission1AtSite[pos + 1],
                      emission0minus1AtSite[pos + 1], emission2minus0AtSite[pos + 1], obsIsZero, obsIsHomMinor);
    }
    if (pos % scalingSkip == 0) {
      m_scalingBuffer = 1.f / getSumOfVector(currentBeta);
      currentBeta = elementWiseMultVectorScalar(currentBeta, m_scalingBuffer);
    }
    fillMatrixColumn(beta, currentBeta, pos);
    lastComputedBeta = currentBeta;
    // update distances
    lastGeneticPos = data.geneticPositions[pos];
    lastPhysicalPos = data.physicalPositions[pos];
  }
  return beta;
}

void HMM::getPreviousBeta(float recDistFromPrevious, vector<float>& lastComputedBeta, vector<float>& BL,
                          vector<float>& BU, vector<float>& currentBeta, vector<float>& emission1AtSite,
                          vector<float>& emission0minus1AtSite, vector<float>& emission2minus0AtSite, float obsIsZero,
                          float obsIsHomMinor)
{
  vector<float> vec = vector<float>(states);
  for (int k = 0; k < states; k++) {
    float currentEmission_k =
        emission1AtSite[k] + emission0minus1AtSite[k] * obsIsZero + emission2minus0AtSite[k] * obsIsHomMinor;
    vec[k] = lastComputedBeta[k] * currentEmission_k;
  }
  // compute below table
  float sum = 0;
  const vector<float>& B = m_decodingQuant.Bvectors.at(recDistFromPrevious);
  for (int k = 1; k < states; k++) {
    sum += B[k - 1] * vec[k - 1];
    BL[k] = sum;
  }
  // compute above table
  const vector<float>& U = m_decodingQuant.Uvectors.at(recDistFromPrevious);
  const vector<float>& RR = m_decodingQuant.rowRatioVectors.at(recDistFromPrevious);
  for (int k = states - 2; k >= 0; k--) {
    BU[k] = vec[k + 1] * U[k] + RR[k] * BU[k + 1];
  }
  // put them together
  const vector<float>& D = m_decodingQuant.Dvectors.at(recDistFromPrevious);
  for (int k = 0; k < states; k++) {
    currentBeta[k] = BL[k] + vec[k] * D[k] + BU[k];
  }
}
