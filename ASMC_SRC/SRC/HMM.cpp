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
#include <chrono>
#include <exception>
#include <iostream>
#include <limits>
#include <utility>

#include <Eigen/Dense>

#include "AvxDefinitions.hpp"
#include "HmmUtils.hpp"
#include "MemoryUtils.hpp"
#include "StringUtils.hpp"
#include "Timer.hpp"

using namespace std;


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

// constructor
HMM::HMM(Data& _data, const DecodingQuantities& _decodingQuant, DecodingParams _decodingParams, bool useBatches,
         int _scalingSkip)
    : data(_data), m_decodingQuant(_decodingQuant), decodingParams(_decodingParams), scalingSkip(_scalingSkip),
      noBatches(!useBatches)
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

  // FastSMC currently specifies batch size as an input
  // todo: enable this option for ASMC as well?
  if (decodingParams.GERMLINE) {
    m_batchSize = decodingParams.batchSize;
  }

  for (int i = 0; i < m_batchSize; i++) {
    fromBatch.push_back(0);
    toBatch.push_back(sequenceLength);
  }

  // get state threshold
  stateThreshold = getStateThreshold();
  // probabilityThreshold = (1./decodingQuant->states)*stateThreshold;

  probabilityThreshold = 0.f;
  for (int i = 0; i < stateThreshold; i++) {
    probabilityThreshold += m_decodingQuant.initialStateProb.at(i);
  }

  if (decodingParams.noConditionalAgeEstimates) {
    ageThreshold = states;
  } else {
    ageThreshold = stateThreshold;
  }

  startBatch = sequenceLength;
  endBatch = 0;

  if (decodingParams.doPerPairPosteriorMean && !decodingParams.GERMLINE) {
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

PairObservations HMM::makePairObs(int_least8_t iHap, unsigned int ind1, int_least8_t jHap, unsigned int ind2)
{
  PairObservations ret;
  ret.iHap = iHap;
  ret.jHap = jHap;
  ret.iInd = ind1;
  ret.jInd = ind2;

  //\todo: ideally all calls to makeBits would be in one place, but GERMLINE calls are in addToBatch and runLastBatch
  if (!decodingParams.GERMLINE || noBatches) {
    makeBits(ret, 0, sequenceLength);
  }

  return ret;
}

void HMM::makeBits(PairObservations& obs, unsigned from, unsigned to)
{
  unsigned iInd = obs.iInd;
  unsigned jInd = obs.jInd;
  obs.obsBits =
      asmc::subsetXorVec(obs.iHap == 1 ? data.individuals[iInd].genotype1 : data.individuals[iInd].genotype2,
                         obs.jHap == 1 ? data.individuals[jInd].genotype1 : data.individuals[jInd].genotype2, from, to);
  obs.homMinorBits =
      asmc::subsetAndVec(obs.iHap == 1 ? data.individuals[iInd].genotype1 : data.individuals[iInd].genotype2,
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
  if (decodingParams.doPerPairPosteriorMean && !decodingParams.GERMLINE) {
    if (foutPosteriorMeanPerPair) {
      foutPosteriorMeanPerPair.close();
    }
    foutPosteriorMeanPerPair.openOrExit(outFileRoot + ".perPairPosteriorMeans.gz");
  }
  if (decodingParams.doPerPairMAP && !decodingParams.GERMLINE) {
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

   auto t0 = std::chrono::high_resolution_clock::now();
  Timer timer;

  resetDecoding();

  const vector<Individual>& individuals = data.individuals;
  uint64 lastPercentage = -1;
  uint64 N = individuals.size();
  uint64 pairs = 0, pairsJob = 0;

  if (decodingParams.GERMLINE) {
    //create IBD output file
    if (!decodingParams.BIN_OUT) {
      gzoutIBD = gzopen( (decodingParams.outFileRoot + "." + std::to_string(jobInd) + "." + std::to_string(jobs) + ".FastSMC.ibd.gz").c_str(), "w" );
    } else {
      gzoutIBD = gzopen( (decodingParams.outFileRoot + "." + std::to_string(jobInd) + "." + std::to_string(jobs) + ".FastSMC.bibd.gz").c_str(), "wb" );
      writeBinaryInfoIntoFile();
    }
    return;
  }

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
              PairObservations observations = makePairObs(jHap, j, iHap, i);
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
      PairObservations observations = makePairObs(1, i, 2, i);
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

  auto t1 = std::chrono::high_resolution_clock::now();
  auto ticksDecodeAll = t1 - t0;
  printf("\nDecoded %" PRIu64 " pairs in %.9f seconds.\n", pairsJob, timer.update_time());

  // print some stats
  asmc::printPctTime("forward", ticksForward / ticksDecodeAll);
  asmc::printPctTime("backward", ticksBackward / ticksDecodeAll);
  asmc::printPctTime("combine", ticksCombine / ticksDecodeAll);
  asmc::printPctTime("sumOverPairs", ticksSumOverPairs / ticksDecodeAll);
  asmc::printPctTime("outputPerPair", ticksOutputPerPair / ticksDecodeAll);
  asmc::printPctTime("other",
                     1 - (ticksForward + ticksBackward + ticksCombine + ticksSumOverPairs + ticksOutputPerPair) /
                             ticksDecodeAll);

  finishDecoding();
}

void HMM::writeBinaryInfoIntoFile()
{
  int mode;
  if (!decodingParams.doPerPairMAP && !decodingParams.doPerPairPosteriorMean) {
    mode = 0;
  } else {
    if (decodingParams.doPerPairPosteriorMean && decodingParams.doPerPairMAP) {
      mode = 1;
    } else if (decodingParams.doPerPairPosteriorMean) {
      mode = 2;
    } else if (decodingParams.doPerPairMAP) {
      mode = 3;
    }
  }
  gzwrite(gzoutIBD, (char*)&mode, sizeof(int));
  gzwrite(gzoutIBD, (char*)&data.chrNumber, sizeof(int));
  unsigned int nbInd = data.individuals.size();
  unsigned int lengthFamid;
  unsigned int lengthIid;
  gzwrite(gzoutIBD, (char*)&nbInd, sizeof(unsigned int));
  for (unsigned int i = 0; i < nbInd; i++) {
    lengthFamid = data.FamIDList[i].size();
    gzwrite(gzoutIBD, (char*)&lengthFamid, sizeof(unsigned int));
    gzwrite(gzoutIBD, data.FamIDList[i].c_str(), lengthFamid);
    lengthIid = data.IIDList[i].size();
    gzwrite(gzoutIBD, (char*)&lengthIid, sizeof(unsigned int));
    gzwrite(gzoutIBD, data.IIDList[i].c_str(), lengthIid);
  }
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
        PairObservations observations = makePairObs(iHap, i, jHap, j);
        if (noBatches) {
          decode(observations);
        } else {
          addToBatch(m_observationsBatch, observations);
        }
      }
    }
  } else {
    // this is the same individual; only decode across chromosomes
    PairObservations observations = makePairObs(1, i, 2, i);
    if (noBatches) {
      decode(observations);
    } else {
      addToBatch(m_observationsBatch, observations);
    }
  }
}

void HMM::decodeFromGERMLINE(const uint indivID1, const uint indivID2, const uint fromPosition, const uint toPosition)
{
  const vector<Individual>& individuals = data.individuals;

  // indivID1 & indivID2 are indices of chromosomes corresponding to individuals indivID1/2 and indivID2/2
  assert(indivID1 / 2 < individuals.size());
  assert(indivID2 / 2 < individuals.size());
  assert(fromPosition < sequenceLength);
  assert(toPosition < sequenceLength);

  Timer timerASMC;

  // ID of individual j must be smaller than ID of individual i
  unsigned int jInd = indivID1 / 2;
  unsigned int iInd = indivID2 / 2;

  PairObservations observation = makePairObs(indivID1 % 2 == 0 ? 1 : 2, jInd, indivID2 % 2 == 0 ? 1 : 2, iInd);

  if (noBatches) {
    decode(observation, fromPosition, toPosition);
  } else {
    nbBatch = cpt % static_cast<unsigned long>(m_batchSize);
    fromBatch[nbBatch] = fromPosition;
    toBatch[nbBatch] = toPosition;
    addToBatch(batchObservations, observation);
    cpt++;
  }
  if (cpt % 10000 == 0) {
    cout << "\rnumber of decoded segments: " << cpt << "\t"
         << "\tdetected segments: " << nbSegmentsDetected << flush;
  }
  timeASMC += timerASMC.update_time();
}

uint HMM::getStateThreshold()
{
  uint result = 0u;
  const vector <float>& disc = m_decodingQuant.discretization;

  while (disc[result] < static_cast<float>(decodingParams.time) && result < m_decodingQuant.states) {
    result++;
  }
  return result;
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

void HMM::finishFromGERMLINE(){
//  timerASMC.update_time();

  if (!noBatches) {
    runLastBatch(batchObservations);
  }

  gzclose(gzoutIBD);

  // print some stats (will remove)
//  timeASMC += timerASMC.update_time();
//  double ticksDecodeAll = ticksForward + ticksBackward + ticksCombine + ticksOutputPerPair;
//  //printf("\n\n*** ASMC decoded %Ld pairs in %.3f seconds. ***\n\n", cpt, timeASMC);
//  printf("\n");
//  printPctTime("forward", ticksForward / ticksDecodeAll);
//  printPctTime("backward", ticksBackward / ticksDecodeAll);
//  printPctTime("combine", ticksCombine / ticksDecodeAll);
//  //printPctTime("sumOverPairs", ticksSumOverPairs / ticksDecodeAll);
//  printPctTime("outputPerPair", ticksOutputPerPair / ticksDecodeAll);
//  cout << flush;
}

// add pair to batch and run if we have enough
void HMM::addToBatch(vector<PairObservations>& obsBatch, const PairObservations& observations)
{
  obsBatch.push_back(observations);
  if (static_cast<int>(obsBatch.size()) == m_batchSize) {

    // taking the maximum 'to' position and the minimum 'from' position in the batch
    startBatch = *std::min_element(fromBatch.begin(), fromBatch.end());
    endBatch = *std::max_element(toBatch.begin(), toBatch.end());

    unsigned int from = asmc::getFromPosition(data.geneticPositions, startBatch);
    unsigned int to = asmc::getToPosition(data.geneticPositions, endBatch);

    if (decodingParams.GERMLINE) {
      for (auto& obs : obsBatch) {
        makeBits(obs, from, to);
      }
    }

    // decodeBatch saves posteriors into m_alphaBuffer [sequenceLength x states x m_batchSize]
    decodeBatch(obsBatch, m_batchSize, m_batchSize, from, to);

    augmentSumOverPairs(obsBatch, m_batchSize, m_batchSize);
    if ((decodingParams.doPerPairMAP || decodingParams.doPerPairPosteriorMean) && !decodingParams.GERMLINE) {
      writePerPairOutput(m_batchSize, m_batchSize, obsBatch);
    }

    // reinitializing batch variables
    startBatch = sequenceLength;
    endBatch = 0;
    obsBatch.clear();
  }
}

// complete with leftover pairs
void HMM::runLastBatch(vector<PairObservations>& obsBatch)
{
  if (obsBatch.empty()) {
    return;
  }

  auto actualBatchSize = obsBatch.size();

  // taking the maximum To position and the minimum From position in the batch
  startBatch = *std::min_element(fromBatch.begin(), fromBatch.begin() + actualBatchSize);
  endBatch = *std::max_element(toBatch.begin(), toBatch.begin() + actualBatchSize);

  unsigned int from = asmc::getFromPosition(data.geneticPositions, startBatch);
  unsigned int to = asmc::getToPosition(data.geneticPositions, endBatch);

  if (decodingParams.GERMLINE) {
    for (auto& obs : obsBatch) {
      makeBits(obs, from, to);
    }
  }

  // fill to size divisible by VECX
  while (obsBatch.size() % VECX != 0) {
    obsBatch.push_back(obsBatch.back());
  }

  auto paddedBatchSize = obsBatch.size();

  // decodeBatch saves posteriors into m_alphaBuffer [sequenceLength x states x paddedBatchSize]
  decodeBatch(obsBatch, actualBatchSize, paddedBatchSize, from, to);
  augmentSumOverPairs(obsBatch, actualBatchSize, paddedBatchSize);

  if ((decodingParams.doPerPairMAP || decodingParams.doPerPairPosteriorMean) && !decodingParams.GERMLINE) {
    writePerPairOutput(actualBatchSize, paddedBatchSize, obsBatch);
  }

  obsBatch.clear();
}

// decode a batch
void HMM::decodeBatch(const vector<PairObservations>& obsBatch, const std::size_t actualBatchSize,
                      const std::size_t paddedBatchSize, const unsigned from, const unsigned to)
{

  int curBatchSize = static_cast<int>(obsBatch.size());

  float* obsIsZeroBatch = ALIGNED_MALLOC_FLOATS(sequenceLength * curBatchSize);
  float* obsIsTwoBatch = ALIGNED_MALLOC_FLOATS(sequenceLength * curBatchSize);

  for (long int pos = from; pos < to; pos++) {
    for (int v = 0; v < curBatchSize; v++) {
      obsIsZeroBatch[pos * curBatchSize + v] = (!obsBatch[v].obsBits[pos - from] ? 1.0f : 0.0f);
      obsIsTwoBatch[pos * curBatchSize + v] = (obsBatch[v].homMinorBits[pos - from] ? 1.0f : 0.0f);
    }
  }

  auto t0 = std::chrono::high_resolution_clock::now();

  // run forward
  forwardBatch(obsIsZeroBatch, obsIsTwoBatch, curBatchSize, from, to);

  auto t1 = std::chrono::high_resolution_clock::now();
  ticksForward += t1 - t0;

  // run backward
  backwardBatch(obsIsZeroBatch, obsIsTwoBatch, curBatchSize, from, to);

  auto t2 = std::chrono::high_resolution_clock::now();
  ticksBackward += t2 - t1;

  // combine (alpha * beta), normalize and store
  float* scale = obsIsZeroBatch; // reuse buffer but rename to be less confusing
  memset(scale, 0, sequenceLength * curBatchSize * sizeof(scale[0]));
#ifdef NO_SSE
  for (long int pos = from; pos < to; pos++) {
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v++) {
        long int ind = (pos * states + k) * curBatchSize + v;
        m_alphaBuffer[ind] *= m_betaBuffer[ind];
        scale[pos * curBatchSize + v] += m_alphaBuffer[ind];
      }
    }
  }
  for (long int pos = from; pos < to; pos++) {
    for (int v = 0; v < curBatchSize; v++) {
      scale[pos * curBatchSize + v] = 1.0f / scale[pos * curBatchSize + v];
    }
  }
  for (long int pos = from; pos < to; pos++) {
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v++) {
        m_alphaBuffer[(pos * states + k) * curBatchSize + v] *= scale[pos * curBatchSize + v];
      }
    }
  }
#else
  for (long int pos = from; pos < to; pos++) {
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v += VECX) {
        long int ind = (pos * states + k) * curBatchSize + v;
        FLOAT prod = MULT(LOAD(&m_alphaBuffer[ind]), LOAD(&m_betaBuffer[ind]));
        STORE(&m_alphaBuffer[ind], prod);
        STORE(&scale[pos * curBatchSize + v], ADD(LOAD(&scale[pos * curBatchSize + v]), prod));
      }
    }
  }
  for (long int pos = from; pos < to; pos++) {
    for (int v = 0; v < curBatchSize; v += VECX) {
      long int ind = pos * curBatchSize + v;
      STORE(&scale[ind], RECIPROCAL(LOAD(&scale[ind])));
    }
  }
  for (long int pos = from; pos < to; pos++) {
    for (int k = 0; k < states; k++) {
      for (int v = 0; v < curBatchSize; v += VECX) {
        long int ind = (pos * states + k) * curBatchSize + v;
        STORE(&m_alphaBuffer[ind], MULT(LOAD(&m_alphaBuffer[ind]), LOAD(&scale[pos * curBatchSize + v])));
      }
    }
  }
#endif

  auto t3 = std::chrono::high_resolution_clock::now();
  ticksCombine += t3 - t2;
  ALIGNED_FREE(obsIsZeroBatch);
  ALIGNED_FREE(obsIsTwoBatch);

  if (decodingParams.GERMLINE) {
    writePerPairOutputFastSMC(actualBatchSize, paddedBatchSize, obsBatch);
  }
}



// forward step
void HMM::forwardBatch(const float* obsIsZeroBatch, const float* obsIsTwoBatch, int curBatchSize, const unsigned from,
                       const unsigned to)
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
      m_alphaBuffer[(states * from + k) * curBatchSize + v] = m_decodingQuant.initialStateProb[k] * firstEmission;
    }
  }

  float* sums = AU; // reuse buffer but rename to be less confusing
  float* currentAlpha = &m_alphaBuffer[states * from * curBatchSize];
  asmc::calculateScalingBatch(currentAlpha, m_scalingBuffer, sums, curBatchSize, states);
  asmc::applyScalingBatch(currentAlpha, m_scalingBuffer, curBatchSize, states);

  // Induction Step:
  float lastGeneticPos = data.geneticPositions[from];
  int lastPhysicalPos = data.physicalPositions[from];

  for (long int pos = from + 1; pos < to; pos++) {
    // get distances and rates
    float recDistFromPrevious = asmc::roundMorgans(data.geneticPositions[pos] - lastGeneticPos, precision, minGenetic);
    float currentRecRate = asmc::roundMorgans(data.recRateAtMarker[pos], precision, minGenetic);
    float* previousAlpha = &m_alphaBuffer[(pos - 1) * states * curBatchSize];
    float* nextAlpha = &m_alphaBuffer[pos * states * curBatchSize];

    if (decodingParams.decodingSequence) {
      int physDistFromPreviousMinusOne =
          asmc::roundPhysical(data.physicalPositions[pos] - lastPhysicalPos - 1, precision);
      float recDistFromPreviousMinusOne =
          asmc::roundMorgans(recDistFromPrevious - currentRecRate, precision, minGenetic);
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

    if (pos % scalingSkip == 0) {
      asmc::calculateScalingBatch(nextAlpha, m_scalingBuffer, sums, curBatchSize, states);
      asmc::applyScalingBatch(nextAlpha, m_scalingBuffer, curBatchSize, states);
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
void HMM::backwardBatch(const float* obsIsZeroBatch, const float* obsIsTwoBatch, int curBatchSize, const unsigned from,
                        const unsigned to)
{

  // fill pos=sequenceLenght-1 in beta
  for (int k = 0; k < states; k++) {
    for (int v = 0; v < curBatchSize; v++) {
      m_betaBuffer[((to - 1) * states + k) * curBatchSize + v] = 1.0f;
    }
  }
  float* sums = ALIGNED_MALLOC_FLOATS(curBatchSize);

  float* currentBeta = &m_betaBuffer[states * (to - 1) * curBatchSize];
  asmc::calculateScalingBatch(currentBeta, m_scalingBuffer, sums, curBatchSize, states);
  asmc::applyScalingBatch(currentBeta, m_scalingBuffer, curBatchSize, states);

  // Induction Step:
  float* BL = ALIGNED_MALLOC_FLOATS(curBatchSize);
  float* BU = ALIGNED_MALLOC_FLOATS(states * curBatchSize);
  float* vec = ALIGNED_MALLOC_FLOATS(states * curBatchSize);

  float lastGeneticPos = data.geneticPositions[to - 1];
  int lastPhysicalPos = data.physicalPositions[to - 1];

  for (long int pos = to - 2; pos >= from; pos--) {
    // get distances and rates
    float recDistFromPrevious = asmc::roundMorgans(lastGeneticPos - data.geneticPositions[pos], precision, minGenetic);
    float currentRecRate = asmc::roundMorgans(data.recRateAtMarker[pos], precision, minGenetic);
    float* currentBeta = &m_betaBuffer[pos * states * curBatchSize];
    float* lastComputedBeta = &m_betaBuffer[(pos + 1) * states * curBatchSize];

    if (decodingParams.decodingSequence) {
      int physDistFromPreviousMinusOne =
          asmc::roundPhysical(lastPhysicalPos - data.physicalPositions[pos] - 1, precision);
      float recDistFromPreviousMinusOne = asmc::roundMorgans(recDistFromPrevious - currentRecRate, precision, minGenetic);
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
      asmc::calculateScalingBatch(currentBeta, m_scalingBuffer, sums, curBatchSize, states);
      asmc::applyScalingBatch(currentBeta, m_scalingBuffer, curBatchSize, states);
    }
    // update distances
    lastGeneticPos = data.geneticPositions[pos];
    lastPhysicalPos = data.physicalPositions[pos];
  }

  ALIGNED_FREE(sums);
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

  auto t0 = std::chrono::high_resolution_clock::now();

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

  auto t1 = std::chrono::high_resolution_clock::now();
  ticksSumOverPairs += t1 - t0;
}

float HMM::getPosteriorMean(const vector<float>& posterior)
{

  float normalization = 1.f / std::accumulate(posterior.begin(), posterior.end(), 0.f);

  float posteriorMean = 0.f;
  for (auto k = 0ul; k < posterior.size(); k++) {
    posteriorMean += normalization * posterior[k] * m_decodingQuant.expectedTimes[k];
  }
  return posteriorMean;
}

float HMM::getMAP(vector<float> posterior)
{
  vector<float> ratioPriorPosterior(posterior.size());
  std::transform(posterior.begin(), posterior.end(), m_decodingQuant.initialStateProb.begin(),
                 ratioPriorPosterior.begin(), std::divides<>());
  const auto MAP_location = std::distance(ratioPriorPosterior.begin(),
                                          std::max_element(ratioPriorPosterior.begin(), ratioPriorPosterior.end()));
  return m_decodingQuant.expectedTimes[MAP_location];
}

// write an IBD segment into output file
void HMM::writePairIBD(const PairObservations& obs, unsigned int posStart, unsigned int posEnd, float prob,
                       vector<float>& posterior, int v, int paddedBatchSize)
{
  nbSegmentsDetected++;
  if (!decodingParams.BIN_OUT) {
    string ind1 = data.FamIDList[obs.iInd] + "\t" + data.IIDList[obs.iInd] + "\t" + to_string(obs.iHap) + "\t";
    string ind2 = data.FamIDList[obs.jInd] + "\t" + data.IIDList[obs.jInd] + "\t" + to_string(obs.jHap) + "\t";
    string segment = std::to_string(data.chrNumber) + "\t" + std::to_string(data.physicalPositions[posStart]) + "\t" +
                     std::to_string(data.physicalPositions[posEnd]) + "\t" +
                     std::to_string(prob / (posEnd - posStart + 1));
    gzwrite(gzoutIBD, (ind1 + ind2 + segment).c_str(), (ind1 + ind2 + segment).size());
    if (decodingParams.doPerPairPosteriorMean) {
      float postMean = getPosteriorMean(posterior);
      string post = "\t" + std::to_string(postMean);
      gzwrite(gzoutIBD, post.c_str(), post.size());
    }
    if (decodingParams.doPerPairMAP) {
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
    pos[0] = data.physicalPositions[posStart];
    pos[1] = data.physicalPositions[posEnd];
    float score = prob / (posEnd - posStart + 1);
    gzwrite(gzoutIBD, (char*)&ind[0], sizeof(int));
    gzwrite(gzoutIBD, &hap[0], sizeof(char));
    gzwrite(gzoutIBD, (char*)&ind[1], sizeof(int));
    gzwrite(gzoutIBD, &hap[1], sizeof(char));
    gzwrite(gzoutIBD, (char*)&pos[0], sizeof(unsigned int));
    gzwrite(gzoutIBD, (char*)&pos[1], sizeof(unsigned int));
    gzwrite(gzoutIBD, (char*)&score, sizeof(float));
    if (decodingParams.doPerPairPosteriorMean) {
      float postMean = getPosteriorMean(posterior);
      gzwrite(gzoutIBD, (char*)&postMean, sizeof(float));
    }
    if (decodingParams.doPerPairMAP) {
      float map = getMAP(posterior);
      gzwrite(gzoutIBD, (char*)&map, sizeof(float));
    }
  }
}

void HMM::writePerPairOutputFastSMC(int actualBatchSize, int paddedBatchSize, const vector<PairObservations>& obsBatch)
{
  for (int v = 0; v < actualBatchSize; v++) {
    bool isIBD = false, isIBD1 = false, isIBD2 = false,
         isIBD3 = false; // true if previous position is IBD, false otherwise
    unsigned int startIBD = 0, startIBD1 = 0, startIBD2 = 0, startIBD3 = 0;
    unsigned int endIBD = 0, endIBD1 = 0, endIBD2 = 0, endIBD3 = 0;
    float posteriorIBD = 0;  // cumulative posterior on an IBD segment
    vector<float> posterior;
    vector<float> sum_posterior_per_state;
    vector<float> prev_sum_posterior_per_state;

    if (decodingParams.doPerPairPosteriorMean || decodingParams.doPerPairMAP) {
      for (uint k = 0; k < ageThreshold; k++) {
        posterior.push_back(0);
        sum_posterior_per_state.push_back(0);
        prev_sum_posterior_per_state.push_back(0);
      }
    }

    if (decodingParams.GERMLINE) {
      // remove these 2 lines if you want the preprocessing step to be less permissive
      // TODO : add a flag for this option
      fromBatch[v] = startBatch;
      toBatch[v] = endBatch;
    }

    for (unsigned int pos = fromBatch[v]; pos < toBatch[v]; pos++) {
      float sum = 0;

      if (decodingParams.doPerPairPosteriorMean || decodingParams.doPerPairMAP) {
        for (uint k = 0; k < ageThreshold; k++) {
          float posterior_pos_state_pair = m_alphaBuffer[(pos * states + k) * paddedBatchSize + v];
          posterior[k] = posterior_pos_state_pair;
          prev_sum_posterior_per_state[k] = sum_posterior_per_state[k];
          sum_posterior_per_state[k] += posterior_pos_state_pair;
          if (k < stateThreshold) {
            sum += posterior_pos_state_pair;
          }
        }
      } else {
        for (uint k = 0; k < stateThreshold; k++) {
          float posterior_pos_state_pair = m_alphaBuffer[(pos * states + k) * paddedBatchSize + v];
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
            writePairIBD(obsBatch[v], startIBD, endIBD, posteriorIBD, prev_sum_posterior_per_state, v, paddedBatchSize);
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
            writePairIBD(obsBatch[v], startIBD, endIBD, posteriorIBD, prev_sum_posterior_per_state, v, paddedBatchSize);
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
            writePairIBD(obsBatch[v], startIBD, endIBD, posteriorIBD, prev_sum_posterior_per_state, v, paddedBatchSize);
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
          writePairIBD(obsBatch[v], startIBD1, endIBD1, posteriorIBD, prev_sum_posterior_per_state, v, paddedBatchSize);
          posteriorIBD = 0;
        } else if (isIBD2) {
          endIBD2 = pos - 1;
          writePairIBD(obsBatch[v], startIBD2, endIBD2, posteriorIBD, prev_sum_posterior_per_state, v, paddedBatchSize);
          posteriorIBD = 0;
        } else if (isIBD3) {
          endIBD3 = pos - 1;
          writePairIBD(obsBatch[v], startIBD3, endIBD3, posteriorIBD, prev_sum_posterior_per_state, v, paddedBatchSize);
          posteriorIBD = 0;
        }
        isIBD = false, isIBD1 = false, isIBD2 = false, isIBD3 = false;
      }
    }
  }
}

// will eventually write binary output instead of gzipped
void HMM::writePerPairOutput(int actualBatchSize, int paddedBatchSize, const vector<PairObservations>& obsBatch)
{

  auto t0 = std::chrono::high_resolution_clock::now();

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

  auto t1 = std::chrono::high_resolution_clock::now();
  ticksOutputPerPair += t1 - t0;
}

// *********************************************************************
// non-batched computations (for debugging and pedagogical reasons only)
// *********************************************************************

vector<vector<float>> HMM::decode(const PairObservations& observations) {
  return decode(observations, 0, sequenceLength);
}

vector<vector<float>> HMM::decode(const PairObservations& observations, const unsigned from, const unsigned to)
{

  auto t0 = std::chrono::high_resolution_clock::now();
  vector<vector<float>> forwardOut = forward(observations, from, to);
  auto t1 = std::chrono::high_resolution_clock::now();
  ticksForward += t1 - t0;

  vector<vector<float>> backwardOut = backward(observations, from, to);
  auto t2 = std::chrono::high_resolution_clock::now();
  ticksBackward += t2 - t1;

  vector<vector<float>> posterior = asmc::elementWiseMultMatrixMatrix(forwardOut, backwardOut);
  posterior = asmc::normalizeMatrixColumns(posterior);

  auto t3 = std::chrono::high_resolution_clock::now();
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
vector<vector<float>> HMM::forward(const PairObservations& observations, const unsigned from, const unsigned to)
{

  vector<vector<float>> alpha(states, vector<float>(sequenceLength));

  uint emissionIndex = observations.obsBits[from] ? 1 : 0;
  // if both samples are carriers, there are two distinguished, otherwise, it's the or
  // (use previous xor). This affects the number of undistinguished for the site
  uint distinguished = observations.homMinorBits[from] ? 2 : emissionIndex;
  uint undistinguished = useCSFSatThisPosition[from] ? data.undistinguishedCounts[from][distinguished] : -1;
  vector<float> emission = getEmission(from, distinguished, undistinguished, emissionIndex);

  vector<float> firstAlpha = asmc::elementWiseMultVectorVector(m_decodingQuant.initialStateProb, emission);

  // cumpute scaling (sum of current alpha vector)
  float m_scalingBuffer = 1.f / asmc::getSumOfVector(firstAlpha);
  // normalize current alpha vector to 1
  firstAlpha = asmc::elementWiseMultVectorScalar(firstAlpha, m_scalingBuffer);

  asmc::fillMatrixColumn(alpha, firstAlpha, from);
  // Induction Step:
  vector<float> nextAlpha(states);
  vector<float> alphaC(states + 1);
  vector<float> previousAlpha = firstAlpha;
  float lastGeneticPos = data.geneticPositions[from];
  int lastPhysicalPos = data.physicalPositions[from];

  for (long int pos = from + 1; pos < to; pos++) {
    auto t0 = std::chrono::high_resolution_clock::now();
    // get distances and rates
    float recDistFromPrevious = asmc::roundMorgans(data.geneticPositions[pos] - lastGeneticPos, precision, minGenetic);
    float currentRecRate = asmc::roundMorgans(data.recRateAtMarker[pos], precision, minGenetic);
    // if both samples are carriers, there are two distinguished, otherwise, it's the
    // or (use previous xor). This affects the number of undistinguished for the site
    float obsIsZero = !observations.obsBits[pos] ? 1.0f : 0.0f;
    float obsIsHomMinor = observations.homMinorBits[pos] ? 1.0f : 0.0f;

    if (decodingParams.decodingSequence) {
      int physDistFromPreviousMinusOne = asmc::roundPhysical(data.physicalPositions[pos] - lastPhysicalPos - 1, precision);
      float recDistFromPreviousMinusOne =
          asmc::roundMorgans(recDistFromPrevious - currentRecRate, precision, minGenetic);
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
    auto t1 = std::chrono::high_resolution_clock::now();
    if (pos % scalingSkip == 0) {
      // compute scaling (sum of current alpha vector)
      m_scalingBuffer = 1.f / asmc::getSumOfVector(nextAlpha);
      // normalize current alpha vector to 1
      nextAlpha = asmc::elementWiseMultVectorScalar(nextAlpha, m_scalingBuffer);
    }
    asmc::fillMatrixColumn(alpha, nextAlpha, pos);
    previousAlpha = nextAlpha;
    auto t2 = std::chrono::high_resolution_clock::now();
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
vector<vector<float>> HMM::backward(const PairObservations& observations, const unsigned from, const unsigned to)
{
  vector<vector<float>> beta(states, vector<float>(sequenceLength));

  vector<float> lastBeta(states);
  for (uint i = 0; i < lastBeta.size(); i++) {
    lastBeta[i] = 1.f;
  }
  // normalize current alpha vector to 1
  float m_scalingBuffer = 1.f / asmc::getSumOfVector(lastBeta);
  lastBeta = asmc::elementWiseMultVectorScalar(lastBeta, m_scalingBuffer);
  asmc::fillMatrixColumn(beta, lastBeta, to - 1);
  // Induction Step:
  vector<float> currentBeta(states);
  vector<float> BL(states);
  vector<float> BU(states);
  vector<float> lastComputedBeta = lastBeta;
  float lastGeneticPos = data.geneticPositions[to - 1];
  int lastPhysicalPos = data.physicalPositions[to - 1];

  for (long int pos = to - 2; pos >= from; pos--) {
    // get distances and rates
    float recDistFromPrevious = asmc::roundMorgans(lastGeneticPos - data.geneticPositions[pos], precision, minGenetic);
    float currentRecRate = asmc::roundMorgans(data.recRateAtMarker[pos], precision, minGenetic);
    // if both samples are carriers, there are two distinguished, otherwise, it's the
    // or (use previous xor). This affects the number of undistinguished for the site
    float obsIsZero = !observations.obsBits[pos + 1] ? 1.0f : 0.0f;
    float obsIsHomMinor = observations.homMinorBits[pos + 1] ? 1.0f : 0.0f;
    if (decodingParams.decodingSequence) {
      int physDistFromPreviousMinusOne =
          asmc::roundPhysical(lastPhysicalPos - data.physicalPositions[pos] - 1, precision);
      float recDistFromPreviousMinusOne =
          asmc::roundMorgans(recDistFromPrevious - currentRecRate, precision, minGenetic);
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
      m_scalingBuffer = 1.f / asmc::getSumOfVector(currentBeta);
      currentBeta = asmc::elementWiseMultVectorScalar(currentBeta, m_scalingBuffer);
    }
    asmc::fillMatrixColumn(beta, currentBeta, pos);
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
