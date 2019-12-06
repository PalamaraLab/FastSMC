#ifndef HMM_H
#define HMM_H

#include <algorithm>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <string>
#include <vector>
#include <zlib.h>

#include <emmintrin.h>
#include <math.h>
#include <pmmintrin.h>
#include <xmmintrin.h>

#include "MemoryUtils.hpp"
#include "Timer.hpp"
#include "Types.hpp"
#include <sys/types.h>

#include <sstream>

using namespace std;

// one of these should be defined in Makefile
// #define NO_SSE
// #define SSE
// #define AVX
// #define AVX512

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
uint64 t1sum = 0, t2sum = 0, t1sumBack = 0;
uint64 tstart = 0, ticksForward = 0, ticksBackward = 0, ticksCombine = 0, ticksSumOverPairs = 0, ticksOutputPerPair = 0;
Timer timerASMC;

// gets genotypes for decoding  (xor --> 1 if het)
vector<bool> xorVec(const vector<bool>& x, const vector<bool>& y, unsigned int from, unsigned int to)
{
  vector<bool> ret(to - from);
  for (unsigned int i = from; i < to; i++) {
    ret[i - from] = x[i] ^ y[i];
  }
  return ret;
}

// computes and of genotype, used to distinguish homozygous minor/derived from homozygous major/ancestral
vector<bool> andVec(const vector<bool>& x, const vector<bool>& y, unsigned int from, unsigned int to)
{
  vector<bool> ret(to - from);
  for (unsigned int i = from; i < to; i++) {
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
      cerr << fileName << " should have \"intervalStart\texpectedCoalescentTime\tintervalEnd\" at each line." << endl;
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

// individual ids and XOR/AND of genotypes
struct PairObservations {
  char iHap, jHap;
  unsigned int iInd, jInd;
  vector<bool> obsBits;
  vector<bool> homMinorBits;
};

// does the linear-time decoding
class HMM
{

  float* alphaBuffer;
  float* betaBuffer;
  float* scalingBuffer;
  float* allZeros;

  // this will be needed for decoding
  Data* data;
  const DecodingQuantities* decodingQuant;
  const DecodingParams* decodingParams;

  string outFileRoot;
  string expectedCoalTimesFile;

  long int sequenceLength;
  uint stateThreshold;
  uint ageThreshold;
  int states;

  int scalingSkip;

  vector<float> expectedCoalTimes;
  vector<bool> useCSFSatThisPosition;

  vector<vector<float>> emission1AtSite;
  vector<vector<float>> emission0minus1AtSite;
  vector<vector<float>> emission2minus0AtSite;

  bool noBatches;
  int currPair = 0;

  float timeASMC = 0;
  int batchSize;

  vector<PairObservations> batchObservations;
  unsigned int startBatch;
  unsigned int endBatch;
  vector<unsigned int> fromBatch;
  vector<unsigned int> toBatch;
  unsigned long int cpt = 0;
  unsigned long int nbSegmentsDetected = 0;
  int nbBatch;
  float probabilityThreshold;

  const int precision = 2;
  const double minGenetic = 1e-10f;

  // FILE
  gzFile gzoutIBD;

public:
  // constructor
  HMM()
  {
    data = 0;
    decodingQuant = 0;
    decodingParams = 0;
  }

  void buildHMM(Data& _data, const DecodingQuantities& _decodingQuant, DecodingParams& _decodingParams,
                float timeInitialization = 0, int _scalingSkip = 1)
  {
    timerASMC.update_time();
    timeASMC += timeInitialization;

    cout << "Instruction set: " << MODE << endl;

    data = &_data;
    decodingQuant = &_decodingQuant;
    decodingParams = &_decodingParams;
    batchSize = decodingParams->batchSize;
    scalingSkip = _scalingSkip;
    noBatches = decodingParams->noBatches;
    sequenceLength = data->sites;
    states = _decodingQuant.states;
    outFileRoot = decodingParams->outFileRoot;
    expectedCoalTimesFile = decodingParams->expectedCoalTimesFile;

    useCSFSatThisPosition = vector<bool>(sequenceLength, false);
    emission1AtSite = vector<vector<float>>(sequenceLength, vector<float>(states));
    emission0minus1AtSite = vector<vector<float>>(sequenceLength, vector<float>(states));
    emission2minus0AtSite = vector<vector<float>>(sequenceLength, vector<float>(states));
    prepareEmissions();

    for (int i = 0; i < batchSize; i++) {
      fromBatch.push_back(0);
      toBatch.push_back(sequenceLength);
    }

    // get state threshold
    stateThreshold = getStateThreshold();
    // probabilityThreshold = (1./decodingQuant->states)*stateThreshold;
    probabilityThreshold = 0;
    for (int i = 0; i < stateThreshold; i++) {
      probabilityThreshold += decodingQuant->initialStateProb[i];
    }

    if (decodingParams->noConditionalAgeEstimates) {
      ageThreshold = states;
    } else {
      ageThreshold = stateThreshold;
    }

    // cout << "prob threshold : " << probabilityThreshold << " " << (1./decodingQuant->states)*stateThreshold << endl;

    startBatch = sequenceLength;
    endBatch = 0;

    timeASMC += timerASMC.update_time();
  }

  // build a pair for two individuals
  PairObservations makePairObs(const Individual& iInd, char iHap, unsigned int ind1, const Individual& jInd, char jHap,
                               unsigned int ind2, int from = 0, int to = 0)
  {
    PairObservations ret;
    ret.iHap = iHap;
    ret.jHap = jHap;
    ret.iInd = ind1;
    ret.jInd = ind2;
    if (!decodingParams->GERMLINE || noBatches) {
      ret.obsBits = xorVec(iHap == '1' ? iInd.genotype1 : iInd.genotype2, jHap == '1' ? jInd.genotype1 : jInd.genotype2,
                           0, sequenceLength);
      ret.homMinorBits = andVec(iHap == '1' ? iInd.genotype1 : iInd.genotype2,
                                jHap == '1' ? jInd.genotype1 : jInd.genotype2, 0, sequenceLength);
    }
    return ret;
  }

  void makeBits(PairObservations& obs, unsigned int from, unsigned int to)
  {
    unsigned int iInd = obs.iInd;
    unsigned int jInd = obs.jInd;
    obs.obsBits =
        xorVec(obs.iHap == '1' ? data->individuals[iInd].genotype1 : data->individuals[iInd].genotype2,
               obs.jHap == '1' ? data->individuals[jInd].genotype1 : data->individuals[jInd].genotype2, from, to);
    obs.homMinorBits =
        andVec(obs.iHap == '1' ? data->individuals[iInd].genotype1 : data->individuals[iInd].genotype2,
               obs.jHap == '1' ? data->individuals[jInd].genotype1 : data->individuals[jInd].genotype2, from, to);
  }

  void prepareEmissions()
  {
    if (decodingParams->skipCSFSdistance < std::numeric_limits<float>::infinity()) {
      useCSFSatThisPosition[0] = true;
      float lastGenCSFSwasUsed = 0.f;
      for (uint pos = 1; pos < sequenceLength; pos++) {
        if (data->geneticPositions[pos] - lastGenCSFSwasUsed >= decodingParams->skipCSFSdistance) {
          // this position is a CSFS position
          useCSFSatThisPosition[pos] = true;
          lastGenCSFSwasUsed = data->geneticPositions[pos];
        }
      }
    }
    for (uint pos = 0; pos < sequenceLength; pos++) {
      if (useCSFSatThisPosition[pos]) {
        int undistAtThisSiteFor0dist = data->undistinguishedCounts[pos][0];
        int undistAtThisSiteFor1dist = data->undistinguishedCounts[pos][1];
        int undistAtThisSiteFor2dist = data->undistinguishedCounts[pos][2];
        if (decodingParams->foldData) {
          // working with folded *data
          for (int k = 0; k < states; k++) {
            if (undistAtThisSiteFor1dist >= 0) {
              emission1AtSite[pos][k] = decodingParams->decodingSequence
                                            ? decodingQuant->foldedCSFSmap[undistAtThisSiteFor1dist][1][k]
                                            : decodingQuant->foldedAscertainedCSFSmap[undistAtThisSiteFor1dist][1][k];
            } else {
              emission1AtSite[pos][k] = 0.f;
            }
            emission0minus1AtSite[pos][k] =
                decodingParams->decodingSequence
                    ? (decodingQuant->foldedCSFSmap[undistAtThisSiteFor0dist][0][k] - emission1AtSite[pos][k])
                    : (decodingQuant->foldedAscertainedCSFSmap[undistAtThisSiteFor0dist][0][k] -
                       emission1AtSite[pos][k]);
            if (undistAtThisSiteFor2dist >= 0) {
              emission2minus0AtSite[pos][k] =
                  decodingParams->decodingSequence
                      ? (decodingQuant->foldedCSFSmap[undistAtThisSiteFor2dist][0][k] -
                         decodingQuant->foldedCSFSmap[undistAtThisSiteFor0dist][0][k])
                      : (decodingQuant->foldedAscertainedCSFSmap[undistAtThisSiteFor2dist][0][k] -
                         decodingQuant->foldedAscertainedCSFSmap[undistAtThisSiteFor0dist][0][k]);
            } else {
              emission2minus0AtSite[pos][k] =
                  decodingParams->decodingSequence
                      ? (0 - decodingQuant->foldedCSFSmap[undistAtThisSiteFor0dist][0][k])
                      : (0 - decodingQuant->foldedAscertainedCSFSmap[undistAtThisSiteFor0dist][0][k]);
            }
          }
        } else {
          // working with unfolded *data
          for (int k = 0; k < states; k++) {
            if (undistAtThisSiteFor1dist >= 0) {
              emission1AtSite[pos][k] = decodingParams->decodingSequence
                                            ? decodingQuant->CSFSmap[undistAtThisSiteFor1dist][1][k]
                                            : decodingQuant->ascertainedCSFSmap[undistAtThisSiteFor1dist][1][k];
            } else {
              emission1AtSite[pos][k] = 0.f;
            }
            float emission0AtThisSiteAndState = 0.f;
            if (undistAtThisSiteFor0dist >= 0) {
              emission0AtThisSiteAndState = decodingParams->decodingSequence
                                                ? decodingQuant->CSFSmap[undistAtThisSiteFor0dist][0][k]
                                                : decodingQuant->ascertainedCSFSmap[undistAtThisSiteFor0dist][0][k];
            }
            emission0minus1AtSite[pos][k] = emission0AtThisSiteAndState - emission1AtSite[pos][k];
            if (undistAtThisSiteFor2dist >= 0) {
              int dist = 2;
              int undist = undistAtThisSiteFor2dist;
              if (undistAtThisSiteFor2dist == data->totalSamplesBound - 2) {
                // for monomorphic derived, fold to CSFS[0][0]
                dist = 0;
                undist = 0;
              }
              emission2minus0AtSite[pos][k] =
                  decodingParams->decodingSequence
                      ? (decodingQuant->CSFSmap[undist][dist][k] - emission0AtThisSiteAndState)
                      : (decodingQuant->ascertainedCSFSmap[undist][dist][k] - emission0AtThisSiteAndState);
            } else {
              emission2minus0AtSite[pos][k] = 0 - emission0AtThisSiteAndState;
            }
          }
        }
      } else {
        // this position is not a CSFS position
        for (int k = 0; k < states; k++) {
          emission1AtSite[pos][k] = decodingParams->decodingSequence ? decodingQuant->classicEmissionTable[1][k]
                                                                     : decodingQuant->compressedEmissionTable[1][k];
          emission0minus1AtSite[pos][k] =
              decodingParams->decodingSequence
                  ? (decodingQuant->classicEmissionTable[0][k] - decodingQuant->classicEmissionTable[1][k])
                  : (decodingQuant->compressedEmissionTable[0][k] - decodingQuant->compressedEmissionTable[1][k]);
          // emission2 = emission0
          emission2minus0AtSite[pos][k] = 0.f;
        }
      }
    }
  }

  // Decodes all pairs. Returns a sum of all decoded posteriors (sequenceLength x states).
  void decodeAll(int jobs, int jobInd)
  {

    uint64 t0 = Timer::rdtsc();

    scalingBuffer = ALIGNED_MALLOC_FLOATS(batchSize);
    alphaBuffer = ALIGNED_MALLOC_FLOATS(sequenceLength * states * batchSize);
    betaBuffer = ALIGNED_MALLOC_FLOATS(sequenceLength * states * batchSize);
    allZeros = ALIGNED_MALLOC_FLOATS(sequenceLength * batchSize);
    memset(allZeros, 0, sequenceLength * batchSize * sizeof(allZeros[0]));

    if (!decodingParams->GERMLINE) {

      uint64 lastPercentage = -1;
      uint64 pairs = 0, pairsJob = 0;
      uint64 totPairs;

      int N = data->individuals.size();
      if (!decodingParams->withinOnly) {
        totPairs = 2 * N * N - N;
      } else {
        totPairs = N;
      }

      cout << "*****ASMC will decode " << totPairs << " pairs of individuals.*****" << endl;

      uint64 pairsStart = totPairs * (jobInd - 1) / jobs;
      uint64 pairsEnd = totPairs * jobInd / jobs;
      uint64 totPairsJob = pairsEnd - pairsStart;

      // create IBD output file
      if (!decodingParams->BIN_OUT) {
        gzoutIBD = gzopen(
            (decodingParams->outFileRoot + "." + std::to_string(jobInd) + "." + std::to_string(jobs) + ".asmc.ibd.gz")
                .c_str(),
            "w");
      } else {
        gzoutIBD = gzopen(
            (decodingParams->outFileRoot + "." + std::to_string(jobInd) + "." + std::to_string(jobs) + ".asmc.bibd.gz")
                .c_str(),
            "wb");
        writeBinaryInfoIntoFile();
      }

      for (uint i = 0; i < data->individuals.size(); i++) {
        if (!decodingParams->withinOnly) {
          for (uint j = 0; j < i; j++) {
            // different individuals; decode 2 haps x 2 haps
            for (int iHap = 1; iHap <= 2; iHap++) {
              for (int jHap = 1; jHap <= 2; jHap++) {
                if (pairsStart <= pairs && pairs < pairsEnd) {
                  char iHapChar = iHap == 1 ? '1' : '2';
                  char jHapChar = jHap == 1 ? '1' : '2';
                  PairObservations observations =
                      makePairObs(data->individuals[j], jHapChar, j, data->individuals[i], iHapChar, i);
                  if (noBatches) {
                    decode(observations);
                  } else {
                    addToBatch(batchObservations, observations);
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
          PairObservations observations = makePairObs(data->individuals[i], '1', i, data->individuals[i], '2', i);
          if (noBatches) {
            decode(observations);
          } else {
            addToBatch(batchObservations, observations);
          }
          pairsJob++;
        }
        pairs++;
        currPair = pairs;
        uint64 percentage = 100 * pairsJob / totPairsJob;
        if (percentage != lastPercentage) {
          cout << "\rdecoding progress: " << percentage << "%"
               << "\t" << pairs << flush;
        }
        lastPercentage = percentage;
      }

      if (!noBatches) {
        runLastBatch(batchObservations);
      }

      // close output file
      gzclose(gzoutIBD);

      // dealloc here
      ALIGNED_FREE(betaBuffer);
      ALIGNED_FREE(alphaBuffer);
      ALIGNED_FREE(scalingBuffer);
      ALIGNED_FREE(allZeros);

      uint64 t1 = Timer::rdtsc();
      double ticksDecodeAll = t1 - t0;
      timeASMC += timerASMC.update_time();

      // print some stats (will remove)
      printf("\n\n*** ASMC decoded %Ld pairs in %.3f seconds. ***\n\n", pairsJob, timeASMC);
      printPctTime("forward", ticksForward / ticksDecodeAll);
      printPctTime("backward", ticksBackward / ticksDecodeAll);
      printPctTime("combine", ticksCombine / ticksDecodeAll);
      // printPctTime("sumOverPairs", ticksSumOverPairs / ticksDecodeAll);
      printPctTime("outputPerPair", ticksOutputPerPair / ticksDecodeAll);
      printPctTime("other", 1 - (ticksForward + ticksBackward + ticksCombine + ticksSumOverPairs + ticksOutputPerPair) /
                                    ticksDecodeAll);
      cout << flush;

    } else {

      // create IBD output file
      if (!decodingParams->BIN_OUT) {
        gzoutIBD = gzopen((decodingParams->outFileRoot + "." + std::to_string(jobInd) + "." + std::to_string(jobs) +
                           ".FastSMC.ibd.gz")
                              .c_str(),
                          "w");
      } else {
        gzoutIBD = gzopen((decodingParams->outFileRoot + "." + std::to_string(jobInd) + "." + std::to_string(jobs) +
                           ".FastSMC.bibd.gz")
                              .c_str(),
                          "wb");
        writeBinaryInfoIntoFile();
      }
    }
  }

  void writeBinaryInfoIntoFile()
  {
    int mode;
    if (!decodingParams->doPerPairMAP && !decodingParams->doPerPairPosteriorMean) {
      mode = 0;
    } else {
      if (decodingParams->doPerPairPosteriorMean && decodingParams->doPerPairMAP) {
        mode = 1;
      } else if (decodingParams->doPerPairPosteriorMean) {
        mode = 2;
      } else if (decodingParams->doPerPairMAP) {
        mode = 3;
      }
    }
    gzwrite(gzoutIBD, (char*)&mode, sizeof(int));
    gzwrite(gzoutIBD, (char*)&data->chrNumber, sizeof(int));
    unsigned int nbInd = data->individuals.size();
    unsigned int lengthFamid;
    unsigned int lengthIid;
    gzwrite(gzoutIBD, (char*)&nbInd, sizeof(unsigned int));
    for (unsigned int i = 0; i < nbInd; i++) {
      lengthFamid = data->FamIDList[i].size();
      gzwrite(gzoutIBD, (char*)&lengthFamid, sizeof(unsigned int));
      gzwrite(gzoutIBD, data->FamIDList[i].c_str(), lengthFamid);
      lengthIid = data->IIDList[i].size();
      gzwrite(gzoutIBD, (char*)&lengthIid, sizeof(unsigned int));
      gzwrite(gzoutIBD, data->IIDList[i].c_str(), lengthIid);
    }
  }

  void decodeFromGERMLINE(unsigned int indivID1, unsigned int indivID2, unsigned int fromPosition,
                          unsigned int toPosition)
  {
    timerASMC.update_time();
    // ID of individual j must be smaller than ID of individual i
    unsigned int jInd = indivID1 / 2;
    unsigned int iInd = indivID2 / 2;

    PairObservations observation = makePairObs(data->individuals[jInd], indivID1 % 2 == 0 ? '1' : '2', jInd,
                                               data->individuals[iInd], indivID2 % 2 == 0 ? '1' : '2', iInd);

    if (noBatches) {
      decodeFastSMC(observation, fromPosition, toPosition);
    } else {
      nbBatch = cpt % batchSize;
      fromBatch[nbBatch] = fromPosition;
      toBatch[nbBatch] = toPosition;
      addToBatchFastSMC(batchObservations, observation);
      cpt++;
    }
    if (cpt % 10000 == 0) {
      cout << "\rnumber of decoded segments: " << cpt << "\t"
           << "\tdetected segments: " << nbSegmentsDetected << flush;
    }
    timeASMC += timerASMC.update_time();
  }

  void finishFromGERMLINE()
  {
    timerASMC.update_time();

    if (!noBatches) {
      runLastBatchFastSMC(batchObservations);
    }

    // dealloc here
    ALIGNED_FREE(betaBuffer);
    ALIGNED_FREE(alphaBuffer);
    ALIGNED_FREE(scalingBuffer);
    ALIGNED_FREE(allZeros);

    // close output file
    gzclose(gzoutIBD);

    // print some stats (will remove)
    timeASMC += timerASMC.update_time();
    double ticksDecodeAll = ticksForward + ticksBackward + ticksCombine + ticksOutputPerPair;
    // printf("\n\n*** ASMC decoded %Ld pairs in %.3f seconds. ***\n\n", cpt, timeASMC);
    printf("\n");
    printPctTime("forward", ticksForward / ticksDecodeAll);
    printPctTime("backward", ticksBackward / ticksDecodeAll);
    printPctTime("combine", ticksCombine / ticksDecodeAll);
    // printPctTime("sumOverPairs", ticksSumOverPairs / ticksDecodeAll);
    printPctTime("outputPerPair", ticksOutputPerPair / ticksDecodeAll);
    cout << flush;
  }

  // private:

  // add pair to batch and run if we have enough
  void addToBatch(vector<PairObservations>& obsBatch, const PairObservations& observations)
  {
    obsBatch.push_back(observations);
    if ((int)obsBatch.size() == batchSize) {
      // decodeBatch saves posteriors into alphaBuffer [sequenceLength x states x batchSize]
      decodeBatch(obsBatch);
      // augmentSumOverPairs(obsBatch, batchSize, batchSize);
      writePerPairOutput(batchSize, batchSize, obsBatch);
      obsBatch.clear();
    }
  }

  // add pair to batch and run if we have enough
  void addToBatchFastSMC(vector<PairObservations>& obsBatch, PairObservations& observations)
  {
    obsBatch.push_back(observations);
    if ((int)obsBatch.size() == batchSize) {

      // taking the maximum To position and the minimum From position in the batch
      startBatch = *std::min_element(fromBatch.begin(), fromBatch.end());
      endBatch = *std::max_element(toBatch.begin(), toBatch.end());

      unsigned int from = getFromPosition(startBatch);
      unsigned int to = getToPosition(endBatch);
      for (int i = 0; i < obsBatch.size(); i++) {
        makeBits(obsBatch[i], from, to);
      }

      // decodeBatch saves posteriors into alphaBuffer [sequenceLength x states x batchSize]
      decodeBatchFastSMC(obsBatch, batchSize, batchSize, from, to);

      // reinitializing batch variables
      startBatch = sequenceLength;
      endBatch = 0;
      obsBatch.clear();
    }
  }

  // complete with leftover pairs
  void runLastBatchFastSMC(vector<PairObservations>& obsBatch)
  {
    if (!obsBatch.empty()) {

      int actualBatchSize = obsBatch.size();

      // taking the maximum To position and the minimum From position in the batch
      startBatch = *std::min_element(fromBatch.begin(), fromBatch.begin() + actualBatchSize);
      endBatch = *std::max_element(toBatch.begin(), toBatch.begin() + actualBatchSize);

      unsigned int from = getFromPosition(startBatch);
      unsigned int to = getToPosition(endBatch);

      for (uint i = 0; i < obsBatch.size(); i++) {
        makeBits(obsBatch[i], from, to);
      }

      while (obsBatch.size() % VECX != 0) { // fill to size divisible by VECX
        obsBatch.push_back(obsBatch.back());
      }

      int paddedBatchSize = obsBatch.size();

      // decodeBatch saves posteriors into alphaBuffer [sequenceLength x states x paddedBatchSize]
      decodeBatchFastSMC(obsBatch, actualBatchSize, paddedBatchSize, from, to);
      obsBatch.clear();
    }
  }

  // complete with leftover pairs
  void runLastBatch(vector<PairObservations>& obsBatch)
  {
    if (obsBatch.empty())
      return;
    int actualBatchSize = obsBatch.size();
    while (obsBatch.size() % VECX != 0) // fill to size divisible by VECX
      obsBatch.push_back(obsBatch.back());
    int paddedBatchSize = obsBatch.size();
    // decodeBatch saves posteriors into alphaBuffer [sequenceLength x states x paddedBatchSize]
    decodeBatch(obsBatch);
    // augmentSumOverPairs(obsBatch, actualBatchSize, paddedBatchSize);
    if (decodingParams->doPerPairMAP || decodingParams->doPerPairPosteriorMean) {
      writePerPairOutput(batchSize, batchSize, obsBatch);
    }
    obsBatch.clear();
  }

  // decode a batch
  void decodeBatch(const vector<PairObservations>& obsBatch)
  {

    int curBatchSize = obsBatch.size();

    float* obsIsZeroBatch = ALIGNED_MALLOC_FLOATS(sequenceLength * curBatchSize);
    float* obsIsTwoBatch = ALIGNED_MALLOC_FLOATS(sequenceLength * curBatchSize);

    for (long int pos = 0; pos < sequenceLength; pos++) {
      for (int v = 0; v < curBatchSize; v++) {
        obsIsZeroBatch[pos * curBatchSize + v] = (!obsBatch[v].obsBits[pos] ? 1.0f : 0.0f);
        obsIsTwoBatch[pos * curBatchSize + v] = (obsBatch[v].homMinorBits[pos] ? 1.0f : 0.0f);
      }
    }

    uint64 t0 = Timer::rdtsc();

    // run forward
    forwardBatch(obsIsZeroBatch, obsIsTwoBatch, curBatchSize);

    uint64 t1 = Timer::rdtsc();
    ticksForward += t1 - t0;

    // run backward
    backwardBatch(obsIsZeroBatch, obsIsTwoBatch, curBatchSize);

    uint64 t2 = Timer::rdtsc();
    ticksBackward += t2 - t1;

    // combine (alpha * beta), normalize and store
    float* scale = obsIsZeroBatch; // reuse buffer but rename to be less confusing
    memset(scale, 0, sequenceLength * curBatchSize * sizeof(scale[0]));
#ifdef NO_SSE
    for (long int pos = 0; pos < sequenceLength; pos++) {
      for (int k = 0; k < states; k++) {
        for (int v = 0; v < curBatchSize; v++) {
          long int ind = (pos * states + k) * curBatchSize + v;
          alphaBuffer[ind] *= betaBuffer[ind];
          scale[pos * curBatchSize + v] += alphaBuffer[ind];
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
          alphaBuffer[(pos * states + k) * curBatchSize + v] *= scale[pos * curBatchSize + v];
        }
      }
    }
#else
    for (long int pos = 0; pos < sequenceLength; pos++) {
      for (int k = 0; k < states; k++) {
        for (int v = 0; v < curBatchSize; v += VECX) {
          long int ind = (pos * states + k) * curBatchSize + v;
          FLOAT prod = MULT(LOAD(&alphaBuffer[ind]), LOAD(&betaBuffer[ind]));
          STORE(&alphaBuffer[ind], prod);
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
          STORE(&alphaBuffer[ind], MULT(LOAD(&alphaBuffer[ind]), LOAD(&scale[pos * curBatchSize + v])));
        }
      }
    }

#endif

    uint64 t3 = Timer::rdtsc();
    ticksCombine += t3 - t2;
    ALIGNED_FREE(obsIsZeroBatch);
    ALIGNED_FREE(obsIsTwoBatch);
  }

  // decode a batch
  void decodeBatchFastSMC(const vector<PairObservations>& obsBatch, int actualBatchSize, int paddedBatchSize,
                          unsigned int from, unsigned int to)
  {

    int curBatchSize = obsBatch.size();
    float* obsIsZeroBatch = ALIGNED_MALLOC_FLOATS(sequenceLength * curBatchSize);
    float* obsIsTwoBatch = ALIGNED_MALLOC_FLOATS(sequenceLength * curBatchSize);
    for (unsigned int pos = from; pos < to; pos++) {
      for (int v = 0; v < curBatchSize; v++) {
        obsIsZeroBatch[pos * curBatchSize + v] = (!obsBatch[v].obsBits[pos - from] ? 1.0f : 0.0f);
        obsIsTwoBatch[pos * curBatchSize + v] = (obsBatch[v].homMinorBits[pos - from] ? 1.0f : 0.0f);
      }
    }

    uint64 t0 = Timer::rdtsc();

    // run forward
    forwardBatchFastSMC(obsIsZeroBatch, obsIsTwoBatch, curBatchSize, from, to);

    uint64 t1 = Timer::rdtsc();
    ticksForward += t1 - t0;
    // cout << "ticksForward : " << ticksForward << endl;

    // run backward
    backwardBatchFastSMC(obsIsZeroBatch, obsIsTwoBatch, curBatchSize, from, to);
    uint64 t2 = Timer::rdtsc();
    ticksBackward += t2 - t1;

    // combine (alpha * beta), normalize and store
    float* scale = obsIsZeroBatch; // reuse buffer but rename to be less confusing
    memset(scale, 0, sequenceLength * curBatchSize * sizeof(scale[0]));
#ifdef NO_SSE
    for (unsigned int pos = from; pos < to; pos++) {
      for (int k = 0; k < states; k++) {
        for (int v = 0; v < curBatchSize; v++) {
          long int ind = (pos * states + k) * curBatchSize + v;
          alphaBuffer[ind] *= betaBuffer[ind];
          scale[pos * curBatchSize + v] += alphaBuffer[ind];
        }
      }
    }
    for (unsigned int pos = from; pos < to; pos++) {
      for (int v = 0; v < curBatchSize; v++) {
        scale[pos * curBatchSize + v] = 1.0f / scale[pos * curBatchSize + v];
      }
    }

    for (unsigned int pos = from; pos < to; pos++) {
      for (int k = 0; k < states; k++) {
        for (int v = 0; v < curBatchSize; v++) {
          alphaBuffer[(pos * states + k) * curBatchSize + v] *= scale[pos * curBatchSize + v];
        }
      }
    }
#else
    for (unsigned int pos = from; pos < to; pos++) {
      for (int k = 0; k < states; k++) {
        for (int v = 0; v < curBatchSize; v += VECX) {
          long int ind = (pos * states + k) * curBatchSize + v;

          FLOAT prod = MULT(LOAD(&alphaBuffer[ind]), LOAD(&betaBuffer[ind]));
          STORE(&alphaBuffer[ind], prod);
          STORE(&scale[pos * curBatchSize + v], ADD(LOAD(&scale[pos * curBatchSize + v]), prod));
        }
      }
    }
    for (unsigned int pos = from; pos < to; pos++) {
      for (int v = 0; v < curBatchSize; v += VECX) {
        long int ind = pos * curBatchSize + v;
        STORE(&scale[ind], RECIPROCAL(LOAD(&scale[ind])));
      }
    }
    for (unsigned int pos = from; pos < to; pos++) {
      for (int k = 0; k < states; k++) {
        for (int v = 0; v < curBatchSize; v += VECX) {
          long int ind = (pos * states + k) * curBatchSize + v;
          STORE(&alphaBuffer[ind], MULT(LOAD(&alphaBuffer[ind]), LOAD(&scale[pos * curBatchSize + v])));
        }
      }
    }

#endif

    uint64 t3 = Timer::rdtsc();
    ticksCombine += t3 - t2;
    ALIGNED_FREE(obsIsZeroBatch);
    ALIGNED_FREE(obsIsTwoBatch);
    writePerPairOutput(actualBatchSize, paddedBatchSize, obsBatch);
  }

  // compute scaling factor for an alpha vector
  void scaleBatch(float* alpha, float* scalings, float* sums, int curBatchSize, int pos)
  {
#ifdef NO_SSE
    // compute scaling (sum of current alpha vector)
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v++) {
        sums[v] += alpha[(0 * states + k) * curBatchSize + v];
      }
    }
    for (int v = 0; v < curBatchSize; v++) {
      scalings[v] = 1.0f / sums[v];
    }
#else
    // compute scaling (sum of current alpha vector)
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v += VECX) {
        STORE(&sums[v], ADD(LOAD(&sums[v]), LOAD(&alpha[(0 * states + k) * curBatchSize + v])));
      }
    }
    for (int v = 0; v < curBatchSize; v++) {
      scalings[v] = 1.0f / sums[v];
    }

#endif
  }

  // apply scaling factor to alpha/beta vector
  void applyScaling(float* vec, float* scalings, int curBatchSize, int pos)
  {
#ifdef NO_SSE
    // normalize current alpha vector to 1
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v++) {
        vec[(0 * states + k) * curBatchSize + v] *= scalings[v];
      }
    }
#else
    // normalize current alpha vector to 1
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v += VECX) {
        STORE(&vec[(0 * states + k) * curBatchSize + v],
              MULT(LOAD(&vec[(0 * states + k) * curBatchSize + v]), LOAD(&scalings[v])));
      }
    }

#endif
  }

  // forward step
  void forwardBatchFastSMC(const float* obsIsZeroBatch, const float* obsIsTwoBatch, int curBatchSize, unsigned int from,
                           unsigned int to)
  {

    assert(curBatchSize % VECX == 0);

    float* alphaC = ALIGNED_MALLOC_FLOATS(states * curBatchSize);
    float* AU = ALIGNED_MALLOC_FLOATS(curBatchSize);

    // fill pos=0 in alpha
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v++) {
        float firstEmission = emission1AtSite[from][k] +
                              emission0minus1AtSite[from][k] * obsIsZeroBatch[from * curBatchSize + v] +
                              emission2minus0AtSite[from][k] * obsIsTwoBatch[from * curBatchSize + v];
        alphaBuffer[(states * from + k) * curBatchSize + v] = decodingQuant->initialStateProb[k] * firstEmission;
      }
    }

    float* sums = AU; // reuse buffer but rename to be less confusing
    memset(sums, 0, curBatchSize * sizeof(sums[0]));
    float* currentAlpha = &alphaBuffer[states * from * curBatchSize];
    scaleBatch(currentAlpha, scalingBuffer, sums, curBatchSize, from);
    applyScaling(currentAlpha, scalingBuffer, curBatchSize, from);
    /*scaleBatch(alphaBuffer, scalingBuffer, sums, curBatchSize, from);
    applyScaling(alphaBuffer, scalingBuffer, curBatchSize, from);*/

    // Induction Step:
    double lastGeneticPos = data->geneticPositions[from];
    unsigned long int lastPhysicalPos = data->physicalPositions[from];

    for (long int pos = from + 1; pos < to; pos++) {
      // get distances and rates
      double recDistFromPrevious = roundMorgans(std::max(minGenetic, data->geneticPositions[pos] - lastGeneticPos));
      double currentRecRate = roundMorgans(data->recRateAtMarker[pos]);
      float* previousAlpha = &alphaBuffer[(pos - 1) * states * curBatchSize];
      float* nextAlpha = &alphaBuffer[pos * states * curBatchSize];

      if (decodingParams->decodingSequence) {
        unsigned long int physDistFromPreviousMinusOne =
            roundPhysical(data->physicalPositions[pos] - lastPhysicalPos - 1);
        double recDistFromPreviousMinusOne = roundMorgans(std::max(minGenetic, recDistFromPrevious - currentRecRate));
        vector<float> homozEmission = decodingQuant->homozygousEmissionMap.at(physDistFromPreviousMinusOne);
        getNextAlphaBatched(recDistFromPreviousMinusOne, alphaC, curBatchSize, previousAlpha, pos, allZeros, allZeros,
                            AU, nextAlpha, homozEmission, homozEmission, homozEmission);
        previousAlpha = nextAlpha;
        getNextAlphaBatched(currentRecRate, alphaC, curBatchSize, previousAlpha, pos, obsIsZeroBatch, obsIsTwoBatch, AU,
                            nextAlpha, emission1AtSite[pos], emission0minus1AtSite[pos], emission2minus0AtSite[pos]);
      } else {
        getNextAlphaBatched(recDistFromPrevious, alphaC, curBatchSize, previousAlpha, pos, obsIsZeroBatch,
                            obsIsTwoBatch, AU, nextAlpha, emission1AtSite[pos], emission0minus1AtSite[pos],
                            emission2minus0AtSite[pos]);
      }
      float* sums = AU; // reuse buffer but rename to be less confusing
      memset(sums, 0, curBatchSize * sizeof(sums[0]));

      if (pos % scalingSkip == 0) {
        scaleBatch(nextAlpha, scalingBuffer, sums, curBatchSize, pos);
        applyScaling(nextAlpha, scalingBuffer, curBatchSize, pos);
      }
      // update distances
      lastGeneticPos = data->geneticPositions[pos];
      lastPhysicalPos = data->physicalPositions[pos];
    }

    ALIGNED_FREE(AU);
    ALIGNED_FREE(alphaC);
  }

  // forward step
  void forwardBatch(const float* obsIsZeroBatch, const float* obsIsTwoBatch, int curBatchSize)
  {

    assert(curBatchSize % VECX == 0);

    float* alphaC = ALIGNED_MALLOC_FLOATS(states * curBatchSize);
    float* AU = ALIGNED_MALLOC_FLOATS(curBatchSize);

    // fill pos=0 in alpha
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v++) {
        float firstEmission = emission1AtSite[0][k] + emission0minus1AtSite[0][k] * obsIsZeroBatch[v] +
                              emission2minus0AtSite[0][k] * obsIsTwoBatch[v];
        alphaBuffer[k * curBatchSize + v] = decodingQuant->initialStateProb[k] * firstEmission;
      }
    }

    float* sums = AU; // reuse buffer but rename to be less confusing
    memset(sums, 0, curBatchSize * sizeof(sums[0]));
    scaleBatch(alphaBuffer, scalingBuffer, sums, curBatchSize, 0);
    applyScaling(alphaBuffer, scalingBuffer, curBatchSize, 0);

    // Induction Step:
    double lastGeneticPos = data->geneticPositions[0];
    unsigned long int lastPhysicalPos = data->physicalPositions[0];

    for (long int pos = 1; pos < sequenceLength; pos++) {
      // get distances and rates
      double recDistFromPrevious = roundMorgans(std::max(minGenetic, data->geneticPositions[pos] - lastGeneticPos));
      double currentRecRate = roundMorgans(data->recRateAtMarker[pos]);
      float* previousAlpha = &alphaBuffer[(pos - 1) * states * curBatchSize];
      float* nextAlpha = &alphaBuffer[pos * states * curBatchSize];
      if (decodingParams->decodingSequence) {
        unsigned long int physDistFromPreviousMinusOne =
            roundPhysical(data->physicalPositions[pos] - lastPhysicalPos - 1);
        double recDistFromPreviousMinusOne = roundMorgans(std::max(minGenetic, recDistFromPrevious - currentRecRate));
        vector<float> homozEmission = decodingQuant->homozygousEmissionMap.at(physDistFromPreviousMinusOne);
        getNextAlphaBatched(recDistFromPreviousMinusOne, alphaC, curBatchSize, previousAlpha, pos, allZeros, allZeros,
                            AU, nextAlpha, homozEmission, homozEmission, homozEmission);
        previousAlpha = nextAlpha;
        getNextAlphaBatched(currentRecRate, alphaC, curBatchSize, previousAlpha, pos, obsIsZeroBatch, obsIsTwoBatch, AU,
                            nextAlpha, emission1AtSite[pos], emission0minus1AtSite[pos], emission2minus0AtSite[pos]);
      } else {
        getNextAlphaBatched(recDistFromPrevious, alphaC, curBatchSize, previousAlpha, pos, obsIsZeroBatch,
                            obsIsTwoBatch, AU, nextAlpha, emission1AtSite[pos], emission0minus1AtSite[pos],
                            emission2minus0AtSite[pos]);
      }
      float* sums = AU; // reuse buffer but rename to be less confusing
      memset(sums, 0, curBatchSize * sizeof(sums[0]));
      if (pos % scalingSkip == 0) {
        scaleBatch(nextAlpha, scalingBuffer, sums, curBatchSize, pos);
        applyScaling(nextAlpha, scalingBuffer, curBatchSize, pos);
      }
      // update distances
      lastGeneticPos = data->geneticPositions[pos];
      lastPhysicalPos = data->physicalPositions[pos];
    }

    ALIGNED_FREE(AU);
    ALIGNED_FREE(alphaC);
  }

  // compute next alpha vector in linear time
  void getNextAlphaBatched(float recDistFromPrevious, float* alphaC, int curBatchSize, const float* previousAlpha,
                           uint pos, const float* obsIsZeroBatch, const float* obsIsTwoBatch, float* AU,
                           float* nextAlpha, const vector<float>& emission1AtSite,
                           const vector<float>& emission0minus1AtSite, const vector<float>& emission2minus0AtSite)
  {

    const float* B = &decodingQuant->Bvectors.at(recDistFromPrevious)[0];
    const float* U = &decodingQuant->Uvectors.at(recDistFromPrevious)[0];
    const float* D = &decodingQuant->Dvectors.at(recDistFromPrevious)[0];

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
          AU[v] = U[k - 1] * previousAlpha[(k - 1) * curBatchSize + v] + decodingQuant->columnRatios[k - 1] * AU[v];
        float term = AU[v] + D[k] * previousAlpha[k * curBatchSize + v];
        if (k < states - 1) {
          term += B[k] * alphaC[(k + 1) * curBatchSize + v];
        }
        float currentEmission_k = emission1AtSite[k] +
                                  emission0minus1AtSite[k] * obsIsZeroBatch[pos * curBatchSize + v] +
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
        colRatios_km1 = LOAD1(decodingQuant->columnRatios[k - 1]);
#else
        Ukm1 = LOAD1(&U[k - 1]);
        colRatios_km1 = LOAD1(&decodingQuant->columnRatios[k - 1]);
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
  void backwardBatchFastSMC(const float* obsIsZeroBatch, const float* obsIsTwoBatch, int curBatchSize,
                            unsigned int from, unsigned int to)
  {

    // fill pos=sequenceLenght-1 in beta
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v++) {
        betaBuffer[((to - 1) * states + k) * curBatchSize + v] = 1.0f;
      }
    }
    float* sums = ALIGNED_MALLOC_FLOATS(curBatchSize);
    memset(sums, 0, curBatchSize * sizeof(sums[0]));

    float* currentBeta = &betaBuffer[states * (to - 1) * curBatchSize];
    scaleBatch(currentBeta, scalingBuffer, sums, curBatchSize, to - 1);
    applyScaling(currentBeta, scalingBuffer, curBatchSize, to - 1);
    /*scaleBatch(betaBuffer, scalingBuffer, sums, curBatchSize, to-1);
    applyScaling(betaBuffer, scalingBuffer, curBatchSize, to-1);*/

    // Induction Step:
    float* BL = ALIGNED_MALLOC_FLOATS(curBatchSize);
    float* BU = ALIGNED_MALLOC_FLOATS(states * curBatchSize);
    float* vec = ALIGNED_MALLOC_FLOATS(states * curBatchSize);

    double lastGeneticPos = data->geneticPositions[to - 1];
    unsigned long int lastPhysicalPos = data->physicalPositions[to - 1];

    for (long int pos = to - 2; pos >= from; pos--) {
      // get distances and rates
      double recDistFromPrevious = roundMorgans(std::max(minGenetic, lastGeneticPos - data->geneticPositions[pos]));
      double currentRecRate = roundMorgans(data->recRateAtMarker[pos]);
      float* currentBeta = &betaBuffer[pos * states * curBatchSize];
      float* lastComputedBeta = &betaBuffer[(pos + 1) * states * curBatchSize];

      if (decodingParams->decodingSequence) {
        unsigned long int physDistFromPreviousMinusOne =
            roundPhysical(lastPhysicalPos - data->physicalPositions[pos] - 1);
        double recDistFromPreviousMinusOne = roundMorgans(std::max(minGenetic, recDistFromPrevious - currentRecRate));
        vector<float> homozEmission = decodingQuant->homozygousEmissionMap.at(physDistFromPreviousMinusOne);
        getPreviousBetaBatched(recDistFromPreviousMinusOne, curBatchSize, lastComputedBeta, pos, allZeros, allZeros,
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
        scaleBatch(currentBeta, scalingBuffer, sums, curBatchSize, pos);
        applyScaling(currentBeta, scalingBuffer, curBatchSize, pos);
      }

      // update distances
      lastGeneticPos = data->geneticPositions[pos];
      lastPhysicalPos = data->physicalPositions[pos];
    }

    ALIGNED_FREE(sums);
    ALIGNED_FREE(vec);
    ALIGNED_FREE(BU);
    ALIGNED_FREE(BL);
  }

  // backward step
  void backwardBatch(const float* obsIsZeroBatch, const float* obsIsTwoBatch, int curBatchSize)
  {

    // fill pos=sequenceLenght-1 in beta
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v++) {
        betaBuffer[((sequenceLength - 1) * states + k) * curBatchSize + v] = 1.0f;
      }
    }
    float* sums = ALIGNED_MALLOC_FLOATS(curBatchSize);
    memset(sums, 0, curBatchSize * sizeof(sums[0]));
    scaleBatch(betaBuffer, scalingBuffer, sums, curBatchSize, sequenceLength - 1);
    applyScaling(betaBuffer, scalingBuffer, curBatchSize, sequenceLength - 1);

    // Induction Step:
    float* BL = ALIGNED_MALLOC_FLOATS(curBatchSize);
    float* BU = ALIGNED_MALLOC_FLOATS(states * curBatchSize);
    float* vec = ALIGNED_MALLOC_FLOATS(states * curBatchSize);

    double lastGeneticPos = data->geneticPositions[sequenceLength - 1];
    unsigned long int lastPhysicalPos = data->physicalPositions[sequenceLength - 1];

    for (long int pos = sequenceLength - 2; pos >= 0; pos--) {
      // get distances and rates
      double recDistFromPrevious = roundMorgans(std::max(minGenetic, lastGeneticPos - data->geneticPositions[pos]));
      double currentRecRate = roundMorgans(data->recRateAtMarker[pos]);
      float* currentBeta = &betaBuffer[pos * states * curBatchSize];
      float* lastComputedBeta = &betaBuffer[(pos + 1) * states * curBatchSize];
      if (decodingParams->decodingSequence) {
        unsigned long int physDistFromPreviousMinusOne =
            roundPhysical(lastPhysicalPos - data->physicalPositions[pos] - 1);
        double recDistFromPreviousMinusOne = roundMorgans(std::max(minGenetic, recDistFromPrevious - currentRecRate));
        vector<float> homozEmission = decodingQuant->homozygousEmissionMap.at(physDistFromPreviousMinusOne);
        getPreviousBetaBatched(recDistFromPreviousMinusOne, curBatchSize, lastComputedBeta, pos, allZeros, allZeros,
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
        scaleBatch(currentBeta, scalingBuffer, sums, curBatchSize, sequenceLength - 1);
        applyScaling(currentBeta, scalingBuffer, curBatchSize, sequenceLength - 1);
      }
      // update distances
      lastGeneticPos = data->geneticPositions[pos];
      lastPhysicalPos = data->physicalPositions[pos];
    }

    ALIGNED_FREE(sums);
    ALIGNED_FREE(vec);
    ALIGNED_FREE(BU);
    ALIGNED_FREE(BL);
  }

  // compute previous beta vector in linear time
  void getPreviousBetaBatched(float recDistFromPrevious, int curBatchSize, const float* lastComputedBeta, int pos,
                              const float* obsIsZeroBatch, const float* obsIsTwoBatch, float* vec, float* BU, float* BL,
                              float* currentBeta, const vector<float>& emission1AtSite,
                              const vector<float>& emission0minus1AtSite, const vector<float>& emission2minus0AtSite)
  {
    const vector<float>& B = decodingQuant->Bvectors.at(recDistFromPrevious);
    const vector<float>& U = decodingQuant->Uvectors.at(recDistFromPrevious);
    const vector<float>& RR = decodingQuant->rowRatioVectors.at(recDistFromPrevious);
    const vector<float>& D = decodingQuant->Dvectors.at(recDistFromPrevious);
#ifdef NO_SSE

    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v++) {
        // cout << k << " " << (pos + 1) << " " << endl;
        float currentEmission_k = emission1AtSite[k] +
                                  emission0minus1AtSite[k] * obsIsZeroBatch[(pos + 1) * curBatchSize + v] +
                                  emission2minus0AtSite[k] * obsIsTwoBatch[(pos + 1) * curBatchSize + v];
        // if (v == 0) cout << " " << currentEmission_k;
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
    for (int k = states - 2; k >= 0; k--) {
      for (int v = 0; v < curBatchSize; v++) {
        BU[k * curBatchSize + v] = U[k] * vec[(k + 1) * curBatchSize + v] + RR[k] * BU[(k + 1) * curBatchSize + v];
      }
    }
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
        if (k) {
          BL[v] += B[k - 1] * vec[(k - 1) * curBatchSize + v];
        }
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

  pair<unsigned int, unsigned int> getMAPinterval(vector<float> posterior)
  {
    vector<float> ratioPriorPosterior(posterior.size());
    std::transform(posterior.begin(), posterior.end(), decodingQuant->initialStateProb.begin(),
                   ratioPriorPosterior.begin(), std::divides<float>());
    int MAP = std::distance(ratioPriorPosterior.begin(),
                            std::max_element(ratioPriorPosterior.begin(), ratioPriorPosterior.end()));
    pair<unsigned int, unsigned int> interval = {decodingQuant->discretization[MAP],
                                                 decodingQuant->discretization[MAP + 1]};
    return interval;
  }

  float getMAP(vector<float> posterior)
  {
    vector<float> ratioPriorPosterior(posterior.size());
    std::transform(posterior.begin(), posterior.end(), decodingQuant->initialStateProb.begin(),
                   ratioPriorPosterior.begin(), std::divides<float>());
    int MAP = std::distance(ratioPriorPosterior.begin(),
                            std::max_element(ratioPriorPosterior.begin(), ratioPriorPosterior.end()));
    return decodingQuant->expectedTimes[MAP];
  }

  float getPosteriorMean(vector<float> posterior)
  {
    float posteriorMean = 0;
    float normalization = std::accumulate(posterior.begin(), posterior.end(), 0.0f);
    for (int k = 0; k < posterior.size(); k++) {
      posteriorMean += (1 / normalization) * posterior[k] * decodingQuant->expectedTimes[k];
    }
    return posteriorMean;
  }

  // write an IBD segment into output file
  void writePairIBD(const PairObservations& obs, unsigned int posStart, unsigned int posEnd, float prob,
                    vector<float>& posterior, int v, int paddedBatchSize)
  {
    nbSegmentsDetected++;
    if (!decodingParams->BIN_OUT) {
      string ind1 = data->FamIDList[obs.iInd] + "\t" + data->IIDList[obs.iInd] + "\t" + obs.iHap + "\t";
      string ind2 = data->FamIDList[obs.jInd] + "\t" + data->IIDList[obs.jInd] + "\t" + obs.jHap + "\t";
      string segment = std::to_string(data->chrNumber) + "\t" + std::to_string(data->physicalPositions[posStart]) +
                       "\t" + std::to_string(data->physicalPositions[posEnd]) + "\t" +
                       std::to_string(prob / (posEnd - posStart + 1));
      gzwrite(gzoutIBD, (ind1 + ind2 + segment).c_str(), (ind1 + ind2 + segment).size());
      if (decodingParams->doPerPairPosteriorMean) {
        float postMean = getPosteriorMean(posterior);
        string post = "\t" + std::to_string(postMean);
        gzwrite(gzoutIBD, post.c_str(), post.size());
      }
      if (decodingParams->doPerPairMAP) {
        float map = getMAP(posterior);
        string map_string = "\t" + std::to_string(map);
        gzwrite(gzoutIBD, map_string.c_str(), map_string.size());
      }
      string end_string = "\n";
      gzwrite(gzoutIBD, end_string.c_str(), end_string.size());

    } else {
      int ind[2];
      ind[0] = obs.iInd;
      ind[1] = obs.jInd;
      char hap[2];
      hap[0] = obs.iHap;
      hap[1] = obs.jHap;
      unsigned int pos[2];
      pos[0] = data->physicalPositions[posStart];
      pos[1] = data->physicalPositions[posEnd];
      float score = prob / (posEnd - posStart + 1);
      gzwrite(gzoutIBD, (char*)&ind[0], sizeof(int));
      gzwrite(gzoutIBD, &hap[0], sizeof(char));
      gzwrite(gzoutIBD, (char*)&ind[1], sizeof(int));
      gzwrite(gzoutIBD, &hap[1], sizeof(char));
      gzwrite(gzoutIBD, (char*)&pos[0], sizeof(unsigned int));
      gzwrite(gzoutIBD, (char*)&pos[1], sizeof(unsigned int));
      gzwrite(gzoutIBD, (char*)&score, sizeof(float));
      if (decodingParams->doPerPairPosteriorMean) {
        float postMean = getPosteriorMean(posterior);
        gzwrite(gzoutIBD, (char*)&postMean, sizeof(float));
      }
      if (decodingParams->doPerPairMAP) {
        float map = getMAP(posterior);
        gzwrite(gzoutIBD, (char*)&map, sizeof(float));
      }
    }
  }

  void writePerPairOutput(int actualBatchSize, int paddedBatchSize, const vector<PairObservations>& obsBatch)
  {

    uint64 t0 = Timer::rdtsc();
    for (int v = 0; v < actualBatchSize; v++) {
      bool isIBD = false, isIBD1 = false, isIBD2 = false,
           isIBD3 = false; // true if previous position is IBD, false otherwise
      unsigned int startIBD = 0, startIBD1 = 0, startIBD2 = 0, startIBD3 = 0;
      unsigned int endIBD = 0, endIBD1 = 0, endIBD2 = 0, endIBD3 = 0;
      float posteriorSite = 0; // posterior on a site
      float posteriorIBD = 0;  // cumulative posterior on an IBD segment
      vector<float> posterior;
      vector<float> sum_posterior_per_state;
      vector<float> prev_sum_posterior_per_state;

      if (decodingParams->doPerPairPosteriorMean || decodingParams->doPerPairMAP) {
        for (uint k = 0; k < ageThreshold; k++) {
          posterior.push_back(0);
          sum_posterior_per_state.push_back(0);
          prev_sum_posterior_per_state.push_back(0);
        }
      }

      if (decodingParams->GERMLINE) {
        // remove these 2 lines if you want the preprocessing step to be less permissive
        // TODO : add a flag for this option
        fromBatch[v] = startBatch;
        toBatch[v] = endBatch;
      }

      for (unsigned int pos = fromBatch[v]; pos < toBatch[v]; pos++) {
        float sum = 0;

        if (decodingParams->doPerPairPosteriorMean || decodingParams->doPerPairMAP) {
          for (uint k = 0; k < ageThreshold; k++) {
            float posterior_pos_state_pair = alphaBuffer[(pos * states + k) * paddedBatchSize + v];
            posterior[k] = posterior_pos_state_pair;
            prev_sum_posterior_per_state[k] = sum_posterior_per_state[k];
            sum_posterior_per_state[k] += posterior_pos_state_pair;
            if (k < stateThreshold) {
              sum += posterior_pos_state_pair;
            }
          }
        } else {
          for (uint k = 0; k < stateThreshold; k++) {
            float posterior_pos_state_pair = alphaBuffer[(pos * states + k) * paddedBatchSize + v];
            sum += posterior_pos_state_pair;
          }
        }

        if (sum >= 1000 * probabilityThreshold) {
          if (!isIBD) {
            startIBD = pos;
            sum_posterior_per_state = posterior;
            if (pos > fromBatch[v] && isIBD1) {
              endIBD1 = pos - 1;
              writePairIBD(obsBatch[v], startIBD1, endIBD1, posteriorIBD, prev_sum_posterior_per_state, v,
                           paddedBatchSize);
            } else if (pos > fromBatch[v] && isIBD2) {
              endIBD2 = pos - 1;
              writePairIBD(obsBatch[v], startIBD2, endIBD2, posteriorIBD, prev_sum_posterior_per_state, v,
                           paddedBatchSize);
            } else if (pos > fromBatch[v] && isIBD3) {
              endIBD3 = pos - 1;
              writePairIBD(obsBatch[v], startIBD3, endIBD3, posteriorIBD, prev_sum_posterior_per_state, v,
                           paddedBatchSize);
            }
            posteriorIBD = sum;
          } else {
            posteriorIBD += sum;
          }
          if (pos == toBatch[v] - 1) {
            endIBD = toBatch[v] - 1;
            writePairIBD(obsBatch[v], startIBD, endIBD, posteriorIBD, sum_posterior_per_state, v, paddedBatchSize);
            posteriorIBD = 0;
          }
          isIBD = true;
          isIBD1 = false, isIBD2 = false, isIBD3 = false;
        } else if (sum >= 100 * probabilityThreshold) {
          if (!isIBD1) {
            startIBD1 = pos;
            sum_posterior_per_state = posterior;
            if (pos > fromBatch[v] && isIBD) {
              endIBD = pos - 1;
              writePairIBD(obsBatch[v], startIBD, endIBD, posteriorIBD, prev_sum_posterior_per_state, v,
                           paddedBatchSize);
            } else if (pos > fromBatch[v] && isIBD2) {
              endIBD2 = pos - 1;
              writePairIBD(obsBatch[v], startIBD2, endIBD2, posteriorIBD, prev_sum_posterior_per_state, v,
                           paddedBatchSize);
            } else if (pos > fromBatch[v] && isIBD3) {
              endIBD3 = pos - 1;
              writePairIBD(obsBatch[v], startIBD3, endIBD3, posteriorIBD, prev_sum_posterior_per_state, v,
                           paddedBatchSize);
            }
            posteriorIBD = sum;
          } else {
            posteriorIBD += sum;
          }
          if (pos == toBatch[v] - 1) {
            endIBD1 = toBatch[v] - 1;
            writePairIBD(obsBatch[v], startIBD1, endIBD1, posteriorIBD, sum_posterior_per_state, v, paddedBatchSize);
            posteriorIBD = 0;
          }
          isIBD = false, isIBD2 = false, isIBD3 = false;
          isIBD1 = true;
        } else if (sum >= 10 * probabilityThreshold) {
          if (!isIBD2) {
            startIBD2 = pos;
            sum_posterior_per_state = posterior;
            if (pos > fromBatch[v] && isIBD1) {
              endIBD1 = pos - 1;
              writePairIBD(obsBatch[v], startIBD1, endIBD1, posteriorIBD, prev_sum_posterior_per_state, v,
                           paddedBatchSize);
            } else if (pos > fromBatch[v] && isIBD) {
              endIBD = pos - 1;
              writePairIBD(obsBatch[v], startIBD, endIBD, posteriorIBD, prev_sum_posterior_per_state, v,
                           paddedBatchSize);
            } else if (pos > fromBatch[v] && isIBD3) {
              endIBD3 = pos - 1;
              writePairIBD(obsBatch[v], startIBD3, endIBD3, posteriorIBD, prev_sum_posterior_per_state, v,
                           paddedBatchSize);
            }
            posteriorIBD = sum;
          } else {
            posteriorIBD += sum;
          }
          if (pos == toBatch[v] - 1) {
            endIBD2 = toBatch[v] - 1;
            writePairIBD(obsBatch[v], startIBD2, endIBD2, posteriorIBD, sum_posterior_per_state, v, paddedBatchSize);
            posteriorIBD = 0;
          }
          isIBD = false, isIBD1 = false, isIBD3 = false;
          isIBD2 = true;
        } else if (sum >= probabilityThreshold) {
          if (!isIBD3) {
            startIBD3 = pos;
            sum_posterior_per_state = posterior;
            if (pos > fromBatch[v] && isIBD1) {
              endIBD1 = pos - 1;
              writePairIBD(obsBatch[v], startIBD1, endIBD1, posteriorIBD, prev_sum_posterior_per_state, v,
                           paddedBatchSize);
            } else if (pos > fromBatch[v] && isIBD) {
              endIBD = pos - 1;
              writePairIBD(obsBatch[v], startIBD, endIBD, posteriorIBD, prev_sum_posterior_per_state, v,
                           paddedBatchSize);
            } else if (pos > fromBatch[v] && isIBD2) {
              endIBD2 = pos - 1;
              writePairIBD(obsBatch[v], startIBD2, endIBD2, posteriorIBD, prev_sum_posterior_per_state, v,
                           paddedBatchSize);
            }
            posteriorIBD = sum;
          } else {
            posteriorIBD += sum;
          }
          if (pos == toBatch[v] - 1) {
            endIBD3 = toBatch[v] - 1;
            writePairIBD(obsBatch[v], startIBD3, endIBD3, posteriorIBD, sum_posterior_per_state, v, paddedBatchSize);
            posteriorIBD = 0;
          }
          isIBD = false, isIBD1 = false, isIBD2 = false;
          isIBD3 = true;
        } else {
          if (isIBD) {
            endIBD = pos - 1;
            writePairIBD(obsBatch[v], startIBD, endIBD, posteriorIBD, prev_sum_posterior_per_state, v, paddedBatchSize);
            posteriorIBD = 0;
          } else if (isIBD1) {
            endIBD1 = pos - 1;
            writePairIBD(obsBatch[v], startIBD1, endIBD1, posteriorIBD, prev_sum_posterior_per_state, v,
                         paddedBatchSize);
            posteriorIBD = 0;
          } else if (isIBD2) {
            endIBD2 = pos - 1;
            writePairIBD(obsBatch[v], startIBD2, endIBD2, posteriorIBD, prev_sum_posterior_per_state, v,
                         paddedBatchSize);
            posteriorIBD = 0;
          } else if (isIBD3) {
            endIBD3 = pos - 1;
            writePairIBD(obsBatch[v], startIBD3, endIBD3, posteriorIBD, prev_sum_posterior_per_state, v,
                         paddedBatchSize);
            posteriorIBD = 0;
          }
          isIBD = false, isIBD1 = false, isIBD2 = false, isIBD3 = false;
        }
      }
    }

    uint64 t1 = Timer::rdtsc();
    ticksOutputPerPair += t1 - t0;
  }

  // *****************************************************
  // non-batched computations (for debugging, will remove)
  // *****************************************************

  float printVector(const vector<float>& vec)
  {
    float sum = 0.f;
    for (uint i = 0; i < vec.size(); i++) {
      cout << vec[i] << "\t";
    }
    cout << endl;
    return sum;
  }

  float getSumOfVector(const vector<float>& vec)
  {
    float sum = 0.f;
    for (uint i = 0; i < vec.size(); i++) {
      sum += vec[i];
    }
    return sum;
  }

  vector<float> elementWiseMultVectorScalar(const vector<float>& vec, float val)
  {
    vector<float> ret(vec.size());
    for (uint i = 0; i < vec.size(); i++) {
      ret[i] = vec[i] * val;
    }
    return ret;
  }

  vector<float> elementWiseMultVectorVector(const vector<float>& vec, const vector<float>& factors)
  {
    vector<float> ret(vec.size());
    for (uint i = 0; i < vec.size(); i++) {
      ret[i] = vec[i] * factors[i];
    }
    return ret;
  }

  vector<vector<float>> elementWiseMultMatrixMatrix(const vector<vector<float>>& matrix1,
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

  vector<vector<float>> normalizeMatrixColumns(const vector<vector<float>>& matrix)
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

  void fillMatrixColumn(vector<vector<float>>& matrix, const vector<float>& vec, long int pos)
  {
    for (uint i = 0; i < vec.size(); i++) {
      matrix[i][pos] = vec[i];
    }
  }

  // return the position of a site 0.5 cM before
  int getFromPosition(unsigned int& from)
  {
    unsigned int pos = from;
    double cumGenDist = 0;
    while (cumGenDist < 0.5 && pos > 0) {
      pos--;
      cumGenDist += (data->geneticPositions[pos + 1] - data->geneticPositions[pos]) * 100;
    }
    return pos;
  }

  // return the position of a site 0.5 cM after
  int getToPosition(unsigned int& to)
  {
    int pos = to;
    double cumGenDist = 0;
    while (cumGenDist < 0.5 && pos < sequenceLength) {
      cumGenDist += (data->geneticPositions[pos] - data->geneticPositions[pos - 1]) * 100;
      pos++;
    }
    return pos;
  }

  // convert generation threshold into state state threshold
  uint getStateThreshold()
  {
    uint stateThreshold = 0;
    while (decodingQuant->discretization[stateThreshold] < decodingParams->time &&
           stateThreshold < decodingQuant->states) {
      stateThreshold++;
    }
    return stateThreshold;
  }

  vector<vector<float>> decode(const PairObservations& observations)
  {

    uint64 t0 = Timer::rdtsc();
    vector<vector<float>> forwardOut = forward(observations);
    uint64 t1 = Timer::rdtsc();
    ticksForward += t1 - t0;

    vector<vector<float>> backwardOut = backward(observations);
    uint64 t2 = Timer::rdtsc();
    ticksBackward += t2 - t1;

    vector<vector<float>> posterior = elementWiseMultMatrixMatrix(forwardOut, backwardOut);
    posterior = normalizeMatrixColumns(posterior);

    uint64 t3 = Timer::rdtsc();
    ticksCombine += t3 - t2;

    bool isIBD = false;
    int startIBD = 0;
    int endIBD = 0;
    float posteriorSite = 0; // posterior on a single site
    float posteriorIBD = 0;  // cumulative posterior on an IBD segment
    float prob = (1. / decodingQuant->states) * stateThreshold;
    string ind1;
    string ind2;
    string segment;

    // MAP version
    for (long int pos = 0; pos < sequenceLength; pos++) {
      float sum = 0;
      for (uint k = 0; k < stateThreshold; k++) {
        float posterior_pos_state_pair = posterior[k][pos];
        sum += posterior_pos_state_pair;
      }

      if (sum >= probabilityThreshold) {
        posteriorSite = sum;
        posteriorIBD += posteriorSite;
        if (!isIBD) {
          startIBD = pos;
          isIBD = true;
        }
        if (pos == sequenceLength - 1) {
          endIBD = sequenceLength - 1;
          isIBD = false;
          ind1 = data->FamIDList[observations.iInd] + "\t" + data->IIDList[observations.iInd] + "\t" +
                 observations.iHap + "\t";
          ind2 = data->FamIDList[observations.jInd] + "\t" + data->IIDList[observations.jInd] + "\t" +
                 observations.jHap + "\t";
          segment = std::to_string(data->chrNumber) + "\t" + std::to_string(data->physicalPositions[startIBD]) + "\t" +
                    std::to_string(data->physicalPositions[endIBD]) + "\t" +
                    std::to_string(posteriorIBD / (endIBD - startIBD + 1)) + "\n";
          gzwrite(gzoutIBD, (ind1 + ind2 + segment).c_str(), (ind1 + ind2 + segment).size());
          posteriorIBD = 0;
        }
      } else {
        if (isIBD) {
          endIBD = pos - 1;
          isIBD = false;
          ind1 = data->FamIDList[observations.iInd] + "\t" + data->IIDList[observations.iInd] + "\t" +
                 observations.iHap + "\t";
          ind2 = data->FamIDList[observations.jInd] + "\t" + data->IIDList[observations.jInd] + "\t" +
                 observations.jHap + "\t";
          segment = std::to_string(data->chrNumber) + "\t" + std::to_string(data->physicalPositions[startIBD]) + "\t" +
                    std::to_string(data->physicalPositions[endIBD]) + "\t" +
                    std::to_string(posteriorIBD / (endIBD - startIBD + 1)) + "\n";
          gzwrite(gzoutIBD, (ind1 + ind2 + segment).c_str(), (ind1 + ind2 + segment).size());
        }
        posteriorIBD = 0;
      }
    }

    return posterior;
  }

  void decodeFastSMC(const PairObservations& observations, unsigned int& from, unsigned int& to)
  {

    uint64 t0 = Timer::rdtsc();
    vector<vector<float>> forwardOut = forwardFastSMC(observations, from, to);
    uint64 t1 = Timer::rdtsc();
    ticksForward += t1 - t0;

    vector<vector<float>> backwardOut = backwardFastSMC(observations, from, to);
    uint64 t2 = Timer::rdtsc();
    ticksBackward += t2 - t1;

    vector<vector<float>> posterior = elementWiseMultMatrixMatrix(forwardOut, backwardOut);
    posterior = normalizeMatrixColumns(posterior);

    uint64 t3 = Timer::rdtsc();
    ticksCombine += t3 - t2;
    bool isIBD = false;
    int startIBD = 0;
    int endIBD = 0;
    float posteriorSite = 0; // posterior on a single site
    float posteriorIBD = 0;  // cumulative posterior on an IBD segment

    string ind1;
    string ind2;
    string segment;

    // sum version
    for (long int pos = from; pos < to; pos++) {
      float sum = 0;
      for (uint k = 0; k < stateThreshold; k++) {
        float posterior_pos_state_pair = posterior[k][pos];
        sum += posterior_pos_state_pair;
      }

      if (sum >= probabilityThreshold) {
        posteriorSite = sum;
        posteriorIBD += posteriorSite;
        if (!isIBD) {
          startIBD = pos;
          isIBD = true;
        }
        if (pos == to - 1) {
          endIBD = to - 1;
          isIBD = false;
          ind1 = data->FamIDList[observations.iInd] + "\t" + data->IIDList[observations.iInd] + "\t" +
                 observations.iHap + "\t";
          ind2 = data->FamIDList[observations.jInd] + "\t" + data->IIDList[observations.jInd] + "\t" +
                 observations.jHap + "\t";
          segment = std::to_string(data->chrNumber) + "\t" + std::to_string(data->physicalPositions[startIBD]) + "\t" +
                    std::to_string(data->physicalPositions[endIBD]) + "\t" +
                    std::to_string(posteriorIBD / (endIBD - startIBD + 1)) + "\n";
          gzwrite(gzoutIBD, (ind1 + ind2 + segment).c_str(), (ind1 + ind2 + segment).size());
          posteriorIBD = 0;
        }
      } else {
        if (isIBD) {
          endIBD = pos - 1;
          isIBD = false;
          ind1 = data->FamIDList[observations.iInd] + "\t" + data->IIDList[observations.iInd] + "\t" +
                 observations.iHap + "\t";
          ind2 = data->FamIDList[observations.jInd] + "\t" + data->IIDList[observations.jInd] + "\t" +
                 observations.jHap + "\t";
          segment = std::to_string(data->chrNumber) + "\t" + std::to_string(data->physicalPositions[startIBD]) + "\t" +
                    std::to_string(data->physicalPositions[endIBD]) + "\t" +
                    std::to_string(posteriorIBD / (endIBD - startIBD + 1)) + "\n";
          gzwrite(gzoutIBD, (ind1 + ind2 + segment).c_str(), (ind1 + ind2 + segment).size());
        }
        posteriorIBD = 0;
      }
    }
  }

  double roundMorgans(double gen)
  {
    float gene1e10 = gen * 1e10f;
    int L10 = std::max(0, (int)floor(log10(gene1e10)) - precision);
    float factor = pow(10, L10);
    double rounded = round(gene1e10 / factor) * factor;
    return std::max(minGenetic, rounded / 1e10f);
  }

  unsigned long int roundPhysical(unsigned long int phys)
  {
    // Since HMM for sequence uses distance-1, it can be -1.
    if (phys < -1) {
      cerr << "ERROR. Int overflow " << phys << endl;
      exit(1);
    }
    // regularize 0 or -1 to 1. Since HMM for sequence uses distance-1, it can be -1.
    phys = std::max((unsigned long int)1, phys);
    unsigned long int L10 = std::max(0, (int)floor(log10(phys)) - precision);
    unsigned long int factor = (int)pow(10, L10);
    unsigned long int rounded = round(phys / (float)factor) * factor;
    return rounded;
  }

  vector<float> getEmission(int pos, int distinguished, int undistinguished, int emissionIndex)
  {
    vector<float> emission;
    if (!useCSFSatThisPosition[pos]) {
      // this position is not a CSFS position, use compressed or classic
      emission = decodingParams->decodingSequence ? decodingQuant->classicEmissionTable[emissionIndex]
                                                  : decodingQuant->compressedEmissionTable[emissionIndex];
    } else {
      // this position is a CSFS position
      if (decodingParams->foldData) {
        emission = decodingParams->decodingSequence
                       ? decodingQuant->foldedCSFSmap[undistinguished][emissionIndex]
                       : decodingQuant->foldedAscertainedCSFSmap[undistinguished][emissionIndex];
      } else {
        emission = decodingParams->decodingSequence ? decodingQuant->CSFSmap[undistinguished][distinguished]
                                                    : decodingQuant->ascertainedCSFSmap[undistinguished][distinguished];
      }
    }
    return emission;
  }

  // forward step
  vector<vector<float>> forward(const PairObservations& observations)
  {

    vector<vector<float>> alpha(states, vector<float>(sequenceLength));

    uint emissionIndex = observations.obsBits[0] ? 1 : 0;
    // if both samples are carriers, there are two distinguished, otherwise, it's the or (use previous xor). This
    // affects the number of undistinguished for the site
    uint distinguished = observations.homMinorBits[0] ? 2 : emissionIndex;
    uint undistinguished = useCSFSatThisPosition[0] ? data->undistinguishedCounts[0][distinguished] : -1;
    vector<float> emission = getEmission(0, distinguished, undistinguished, emissionIndex);

    vector<float> firstAlpha = elementWiseMultVectorVector(decodingQuant->initialStateProb, emission);

    // cumpute scaling (sum of current alpha vector)
    float scalingBuffer = 1.f / getSumOfVector(firstAlpha);
    // normalize current alpha vector to 1
    firstAlpha = elementWiseMultVectorScalar(firstAlpha, scalingBuffer);

    fillMatrixColumn(alpha, firstAlpha, 0);
    // Induction Step:
    vector<float> nextAlpha(states);
    vector<float> alphaC(states + 1);
    vector<float> previousAlpha = firstAlpha;
    double lastGeneticPos = data->geneticPositions[0];
    unsigned long int lastPhysicalPos = data->physicalPositions[0];

    for (long int pos = 1; pos < sequenceLength; pos++) {
      unsigned long long t0 = Timer::rdtsc();
      // get distances and rates
      double recDistFromPrevious = roundMorgans(std::max(minGenetic, data->geneticPositions[pos] - lastGeneticPos));
      double currentRecRate = roundMorgans(data->recRateAtMarker[pos]);
      // if both samples are carriers, there are two distinguished, otherwise, it's the or (use previous xor). This
      // affects the number of undistinguished for the site
      float obsIsZero = !observations.obsBits[pos] ? 1.0f : 0.0f;
      float obsIsHomMinor = observations.homMinorBits[pos] ? 1.0f : 0.0f;

      // int distinguished = (1 - obsIsZero) + 2 * obsIsHomMinor;

      // cout << (1 - obsIsZero) << " " << obsIsHomMinor << endl;
      if (decodingParams->decodingSequence) {
        unsigned long int physDistFromPreviousMinusOne =
            roundPhysical(data->physicalPositions[pos] - lastPhysicalPos - 1);
        double recDistFromPreviousMinusOne = roundMorgans(std::max(minGenetic, recDistFromPrevious - currentRecRate));
        vector<float> homozEmission = decodingQuant->homozygousEmissionMap.at(physDistFromPreviousMinusOne);
        getNextAlpha(recDistFromPreviousMinusOne, alphaC, previousAlpha, nextAlpha, homozEmission, homozEmission,
                     homozEmission, 0.0f, 0.0f);
        previousAlpha = nextAlpha;
        getNextAlpha(currentRecRate, alphaC, previousAlpha, nextAlpha, emission1AtSite[pos], emission0minus1AtSite[pos],
                     emission2minus0AtSite[pos], obsIsZero, obsIsHomMinor);
      } else {
        getNextAlpha(recDistFromPrevious, alphaC, previousAlpha, nextAlpha, emission1AtSite[pos],
                     emission0minus1AtSite[pos], emission2minus0AtSite[pos], obsIsZero, obsIsHomMinor);
      }
      // printVector(nextAlpha);
      unsigned long long t1 = Timer::rdtsc();
      if (pos % scalingSkip == 0) {
        // compute scaling (sum of current alpha vector)
        scalingBuffer = 1.f / getSumOfVector(nextAlpha);
        // normalize current alpha vector to 1
        nextAlpha = elementWiseMultVectorScalar(nextAlpha, scalingBuffer);
      }
      fillMatrixColumn(alpha, nextAlpha, pos);
      previousAlpha = nextAlpha;
      unsigned long long t2 = Timer::rdtsc();
      t1sum += t1 - t0;
      t2sum += t2 - t1;
      // update distances
      lastGeneticPos = data->geneticPositions[pos];
      lastPhysicalPos = data->physicalPositions[pos];
    }
    return alpha;
  }

  // forward step
  vector<vector<float>> forwardFastSMC(const PairObservations& observations, unsigned int from, unsigned int to)
  {

    // int segmentLength = to-from;
    vector<vector<float>> alpha(states, vector<float>(sequenceLength));

    uint emissionIndex = observations.obsBits[from] ? 1 : 0;
    // if both samples are carriers, there are two distinguished, otherwise, it's the or (use previous xor). This
    // affects the number of undistinguished for the site
    uint distinguished = observations.homMinorBits[from] ? 2 : emissionIndex;
    uint undistinguished = useCSFSatThisPosition[from] ? data->undistinguishedCounts[from][distinguished] : -1;
    vector<float> emission = getEmission(from, distinguished, undistinguished, emissionIndex);

    vector<float> firstAlpha = elementWiseMultVectorVector(decodingQuant->initialStateProb, emission);

    // cumpute scaling (sum of current alpha vector)
    float scalingBuffer = 1.f / getSumOfVector(firstAlpha);
    // normalize current alpha vector to 1
    firstAlpha = elementWiseMultVectorScalar(firstAlpha, scalingBuffer);

    fillMatrixColumn(alpha, firstAlpha, from);
    // Induction Step:
    vector<float> nextAlpha(states);
    vector<float> alphaC(states + 1);
    vector<float> previousAlpha = firstAlpha;
    double lastGeneticPos = data->geneticPositions[from];
    unsigned long int lastPhysicalPos = data->physicalPositions[from];

    for (long int pos = from + 1; pos < to; pos++) {
      unsigned long long t0 = Timer::rdtsc();
      // get distances and rates
      double recDistFromPrevious = roundMorgans(std::max(minGenetic, data->geneticPositions[pos] - lastGeneticPos));
      double currentRecRate = roundMorgans(data->recRateAtMarker[pos]);
      // if both samples are carriers, there are two distinguished, otherwise, it's the or (use previous xor). This
      // affects the number of undistinguished for the site
      float obsIsZero = !observations.obsBits[pos] ? 1.0f : 0.0f;
      float obsIsHomMinor = observations.homMinorBits[pos] ? 1.0f : 0.0f;

      // int distinguished = (1 - obsIsZero) + 2 * obsIsHomMinor;

      // cout << (1 - obsIsZero) << " " << obsIsHomMinor << endl;
      if (decodingParams->decodingSequence) {
        unsigned long int physDistFromPreviousMinusOne =
            roundPhysical(data->physicalPositions[pos] - lastPhysicalPos - 1);
        double recDistFromPreviousMinusOne = roundMorgans(std::max(minGenetic, recDistFromPrevious - currentRecRate));
        vector<float> homozEmission = decodingQuant->homozygousEmissionMap.at(physDistFromPreviousMinusOne);
        getNextAlpha(recDistFromPreviousMinusOne, alphaC, previousAlpha, nextAlpha, homozEmission, homozEmission,
                     homozEmission, 0.0f, 0.0f);
        previousAlpha = nextAlpha;
        getNextAlpha(currentRecRate, alphaC, previousAlpha, nextAlpha, emission1AtSite[pos], emission0minus1AtSite[pos],
                     emission2minus0AtSite[pos], obsIsZero, obsIsHomMinor);
      } else {
        getNextAlpha(recDistFromPrevious, alphaC, previousAlpha, nextAlpha, emission1AtSite[pos],
                     emission0minus1AtSite[pos], emission2minus0AtSite[pos], obsIsZero, obsIsHomMinor);
      }
      // printVector(nextAlpha);
      unsigned long long t1 = Timer::rdtsc();
      if (pos % scalingSkip == 0) {
        // compute scaling (sum of current alpha vector)
        scalingBuffer = 1.f / getSumOfVector(nextAlpha);
        // normalize current alpha vector to 1
        nextAlpha = elementWiseMultVectorScalar(nextAlpha, scalingBuffer);
      }
      fillMatrixColumn(alpha, nextAlpha, pos);
      previousAlpha = nextAlpha;
      unsigned long long t2 = Timer::rdtsc();
      t1sum += t1 - t0;
      t2sum += t2 - t1;
      // update distances
      lastGeneticPos = data->geneticPositions[pos];
      lastPhysicalPos = data->physicalPositions[pos];
    }
    return alpha;
  }

  void getNextAlpha(float recDistFromPrevious, vector<float>& alphaC, vector<float>& previousAlpha,
                    vector<float>& nextAlpha, vector<float>& emission1AtSite, vector<float>& emission0minus1AtSite,
                    vector<float>& emission2minus0AtSite, float obsIsZero, float obsIsHomMinor)
  {
    alphaC[decodingQuant->states - 1] = previousAlpha[decodingQuant->states - 1];
    for (int k = decodingQuant->states - 2; k >= 0; k--) {
      alphaC[k] = alphaC[k + 1] + previousAlpha[k];
    }
    const float* B = &decodingQuant->Bvectors.at(recDistFromPrevious)[0];
    const float* U = &decodingQuant->Uvectors.at(recDistFromPrevious)[0];
    const float* D = &decodingQuant->Dvectors.at(recDistFromPrevious)[0];
    float AUc = 0;
    for (uint k = 0; k < decodingQuant->states; k++) {
      if (k)
        AUc = U[k - 1] * previousAlpha[k - 1] + decodingQuant->columnRatios[k - 1] * AUc;
      float term = AUc + D[k] * previousAlpha[k];
      if (k < decodingQuant->states - 1) // TODO: just extend B and alphaC?
        term += B[k] * alphaC[k + 1];
      float currentEmission_k =
          emission1AtSite[k] + emission0minus1AtSite[k] * obsIsZero + emission2minus0AtSite[k] * obsIsHomMinor;
      nextAlpha[k] = currentEmission_k * term;
    }
  }

  // backward step
  vector<vector<float>> backward(const PairObservations& observations)
  {

    vector<vector<float>> beta(states, vector<float>(sequenceLength));

    vector<float> lastBeta(states);
    for (uint i = 0; i < lastBeta.size(); i++) {
      lastBeta[i] = 1.f;
    }
    // normalize current alpha vector to 1
    float scalingBuffer = 1.f / getSumOfVector(lastBeta);
    lastBeta = elementWiseMultVectorScalar(lastBeta, scalingBuffer);
    fillMatrixColumn(beta, lastBeta, sequenceLength - 1);
    // Induction Step:
    vector<float> currentBeta(states);
    vector<float> BL(states);
    vector<float> BU(states);
    vector<float> lastComputedBeta = lastBeta;
    double lastGeneticPos = data->geneticPositions[sequenceLength - 1];
    unsigned long int lastPhysicalPos = data->physicalPositions[sequenceLength - 1];
    for (long int pos = sequenceLength - 2; pos >= 0; pos--) {
      // get distances and rates
      double recDistFromPrevious = roundMorgans(std::max(minGenetic, lastGeneticPos - data->geneticPositions[pos]));
      double currentRecRate = roundMorgans(data->recRateAtMarker[pos]);
      // if both samples are carriers, there are two distinguished, otherwise, it's the or (use previous xor). This
      // affects the number of undistinguished for the site
      float obsIsZero = !observations.obsBits[pos + 1] ? 1.0f : 0.0f;         // TODO: vectorize
      float obsIsHomMinor = observations.homMinorBits[pos + 1] ? 1.0f : 0.0f; // TODO: vectorize
      if (decodingParams->decodingSequence) {
        unsigned long int physDistFromPreviousMinusOne =
            roundPhysical(lastPhysicalPos - data->physicalPositions[pos] - 1);
        double recDistFromPreviousMinusOne = roundMorgans(std::max(minGenetic, recDistFromPrevious - currentRecRate));
        vector<float> homozEmission = decodingQuant->homozygousEmissionMap.at(physDistFromPreviousMinusOne);
        getPreviousBeta(recDistFromPreviousMinusOne, lastComputedBeta, BL, BU, currentBeta, homozEmission,
                        homozEmission, homozEmission, 0.0f, 0.0f);
        lastComputedBeta = currentBeta;
        getPreviousBeta(currentRecRate, lastComputedBeta, BL, BU, currentBeta, emission1AtSite[pos + 1],
                        emission0minus1AtSite[pos + 1], emission2minus0AtSite[pos + 1], obsIsZero, obsIsHomMinor);
      } else {
        getPreviousBeta(recDistFromPrevious, lastComputedBeta, BL, BU, currentBeta, emission1AtSite[pos + 1],
                        emission0minus1AtSite[pos + 1], emission2minus0AtSite[pos + 1], obsIsZero, obsIsHomMinor);
      }
      // printVector(currentBeta);
      if (pos % scalingSkip == 0) {
        scalingBuffer = 1.f / getSumOfVector(currentBeta);
        currentBeta = elementWiseMultVectorScalar(currentBeta, scalingBuffer);
      }
      fillMatrixColumn(beta, currentBeta, pos);
      lastComputedBeta = currentBeta;
      // update distances
      lastGeneticPos = data->geneticPositions[pos];
      lastPhysicalPos = data->physicalPositions[pos];
    }
    return beta;
  }

  // backward step
  vector<vector<float>> backwardFastSMC(const PairObservations& observations, unsigned int from, unsigned int to)
  {

    // int segmentLength = to-from;
    vector<vector<float>> beta(states, vector<float>(sequenceLength));

    vector<float> lastBeta(states);
    for (uint i = 0; i < lastBeta.size(); i++) {
      lastBeta[i] = 1.f;
    }
    // normalize current alpha vector to 1
    float scalingBuffer = 1.f / getSumOfVector(lastBeta);
    lastBeta = elementWiseMultVectorScalar(lastBeta, scalingBuffer);
    fillMatrixColumn(beta, lastBeta, to - 1);
    // Induction Step:
    vector<float> currentBeta(states);
    vector<float> BL(states);
    vector<float> BU(states);
    vector<float> lastComputedBeta = lastBeta;
    double lastGeneticPos = data->geneticPositions[to - 1];
    unsigned long int lastPhysicalPos = data->physicalPositions[to - 1];

    for (long int pos = to - 2; pos >= from; pos--) {
      // get distances and rates
      double recDistFromPrevious = roundMorgans(std::max(minGenetic, lastGeneticPos - data->geneticPositions[pos]));
      double currentRecRate = roundMorgans(data->recRateAtMarker[pos]);
      // if both samples are carriers, there are two distinguished, otherwise, it's the or (use previous xor). This
      // affects the number of undistinguished for the site
      float obsIsZero = !observations.obsBits[pos + 1] ? 1.0f : 0.0f;         // TODO: vectorize
      float obsIsHomMinor = observations.homMinorBits[pos + 1] ? 1.0f : 0.0f; // TODO: vectorize
      if (decodingParams->decodingSequence) {
        unsigned long int physDistFromPreviousMinusOne =
            roundPhysical(lastPhysicalPos - data->physicalPositions[pos] - 1);
        double recDistFromPreviousMinusOne = roundMorgans(std::max(minGenetic, recDistFromPrevious - currentRecRate));
        vector<float> homozEmission = decodingQuant->homozygousEmissionMap.at(physDistFromPreviousMinusOne);
        getPreviousBeta(recDistFromPreviousMinusOne, lastComputedBeta, BL, BU, currentBeta, homozEmission,
                        homozEmission, homozEmission, 0.0f, 0.0f);
        lastComputedBeta = currentBeta;
        getPreviousBeta(currentRecRate, lastComputedBeta, BL, BU, currentBeta, emission1AtSite[pos + 1],
                        emission0minus1AtSite[pos + 1], emission2minus0AtSite[pos + 1], obsIsZero, obsIsHomMinor);
      } else {
        getPreviousBeta(recDistFromPrevious, lastComputedBeta, BL, BU, currentBeta, emission1AtSite[pos + 1],
                        emission0minus1AtSite[pos + 1], emission2minus0AtSite[pos + 1], obsIsZero, obsIsHomMinor);
      }
      // printVector(currentBeta);
      if (pos % scalingSkip == 0) {
        scalingBuffer = 1.f / getSumOfVector(currentBeta);
        currentBeta = elementWiseMultVectorScalar(currentBeta, scalingBuffer);
      }
      fillMatrixColumn(beta, currentBeta, pos);
      lastComputedBeta = currentBeta;
      // update distances
      lastGeneticPos = data->geneticPositions[pos];
      lastPhysicalPos = data->physicalPositions[pos];
    }
    return beta;
  }

  void getPreviousBeta(float recDistFromPrevious, vector<float>& lastComputedBeta, vector<float>& BL, vector<float>& BU,
                       vector<float>& currentBeta, vector<float>& emission1AtSite, vector<float>& emission0minus1AtSite,
                       vector<float>& emission2minus0AtSite, float obsIsZero, float obsIsHomMinor)
  {
    vector<float> vec = vector<float>(states);
    for (int k = 0; k < states; k++) {
      float currentEmission_k =
          emission1AtSite[k] + emission0minus1AtSite[k] * obsIsZero + emission2minus0AtSite[k] * obsIsHomMinor;
      vec[k] = lastComputedBeta[k] * currentEmission_k;
    }
    // compute below table
    float sum = 0;
    const vector<float>& B = decodingQuant->Bvectors.at(recDistFromPrevious);
    for (int k = 1; k < states; k++) {
      sum += B[k - 1] * vec[k - 1];
      BL[k] = sum;
    }
    // compute above table
    const vector<float>& U = decodingQuant->Uvectors.at(recDistFromPrevious);
    const vector<float>& RR = decodingQuant->rowRatioVectors.at(recDistFromPrevious);
    for (int k = states - 2; k >= 0; k--) {
      BU[k] = vec[k + 1] * U[k] + RR[k] * BU[k + 1];
    }
    // put them together
    const vector<float>& D = decodingQuant->Dvectors.at(recDistFromPrevious);
    for (int k = 0; k < states; k++) {
      currentBeta[k] = BL[k] + vec[k] * D[k] + BU[k];
    }
  }
};
#endif
