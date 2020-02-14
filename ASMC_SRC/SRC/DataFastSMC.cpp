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

#include <iostream>
#include <string>
#include <vector>

#include "FileUtils.hpp"
#include "StringUtils.hpp"
#include "Types.hpp"

#include "DataFastSMC.hpp"
#include <boost/math/distributions/hypergeometric.hpp>

using namespace std;

DataFastSMC::DataFastSMC(string inFileRoot, unsigned int _sites, int _totalSamplesBound, bool foldToMinorAlleles,
                         bool _decodingUsesCSFS, int jobID, int jobs,
                         vector<pair<unsigned long int, double>>& genetic_map)
    : sites(_sites), totalSamplesBound(_totalSamplesBound), decodingUsesCSFS(_decodingUsesCSFS)
{
  sampleSize = countSamplesLines(inFileRoot);

  siteWasFlippedDuringFolding = vector<bool>(sites, false);

  windowSize = ceil(sqrt((2 * pow(sampleSize, 2) - sampleSize) * 2 / jobs)); // length of a square, in term of #ind
  if (windowSize % 2 != 0) {
    windowSize++;
  }
  w_i = 1; // window for individual i
  int cpt_job = 1;
  int cpt_tot_job = 1;
  while (cpt_tot_job < jobID) {
    w_i++;
    cpt_job = cpt_job + 2;
    cpt_tot_job = cpt_tot_job + cpt_job;
  }
  w_j = ceil((float)(cpt_job - (cpt_tot_job - jobID)) / 2); // window for individual j
  is_j_above_diag = (cpt_job - (cpt_tot_job - jobID)) % 2 == 1 ? true : false;

  readSamplesList(inFileRoot, jobID, jobs); // now we have the sample sizes

  // and the number of sites
  for (uint i = 0; i < FamIDList.size(); i++) {
    individuals.push_back(Individual(sites));
  }

  // now read all the data
  readHaps(inFileRoot, foldToMinorAlleles, jobID, jobs, genetic_map);
  // readMap(inFileRoot);
  // may read freq file here
  // make undistinguished counts
  if (decodingUsesCSFS) {
    makeUndistinguished(foldToMinorAlleles);
  }
}

// bruteforce sampling of hypergeometric
int DataFastSMC::sampleHypergeometric(int populationSize, int numberOfSuccesses, unsigned int sampleSize)
{
  // cout << populationSize << " " << numberOfSuccesses << " " << sampleSize << endl;
  if (numberOfSuccesses < 0 || numberOfSuccesses > populationSize) {
    return -1;
  }
  vector<unsigned short> samplingVector;
  samplingVector = vector<unsigned short>(populationSize, 0);
  for (int i = 0; i < numberOfSuccesses; i++) {
    samplingVector[i] = 1;
  }
  std::random_shuffle(samplingVector.begin(), samplingVector.end());
  int ret = 0;
  for (unsigned int i = 0; i < sampleSize; i++) {
    ret += samplingVector[i];
  }
  return ret;
}

void DataFastSMC::makeUndistinguished(bool foldToMinorAlleles)
{
  cout << "Building undistinguished counts...";
  undistinguishedCounts = vector<vector<int>>(derivedAlleleCounts.size());
  for (uint i = 0; i < derivedAlleleCounts.size(); i++) {
    undistinguishedCounts[i] = vector<int>(3);
    unsigned int derivedAlleles = derivedAlleleCounts[i];
    unsigned int totalSamples = totalSamplesCount[i];
    if (decodingUsesCSFS && totalSamplesBound > totalSamples) {
      cerr << "ERROR. SNP number " << i << " has " << totalSamples << " non-missing individuals, but the CSFS requires "
           << totalSamplesBound << endl;
      exit(1);
    }
    int ancestralAlleles = totalSamples - derivedAlleles;
    if (foldToMinorAlleles && derivedAlleles > ancestralAlleles) {
      cerr << "Minor alleles has frequency > 50%. Data is supposed to be folded.\n";
    }
    for (int distinguished = 0; distinguished < 3; distinguished++) {
      // hypergeometric with (derivedAlleles - distinguished) derived alleles, (samples - 2) samples
      // cout << "At " << i << ", hypergeometric with " << (totalSamples - 2) << " samples, " << (derivedAlleles -
      // distinguished) << " derived remaining after removing 2 samples with " << distinguished << " derived. We have "
      // << (totalSamplesBound - 2) << " undistinguished to draw" << endl;
      int undist = sampleHypergeometric(totalSamples - 2, derivedAlleles - distinguished, totalSamplesBound - 2);
      if (foldToMinorAlleles && (undist + distinguished > totalSamplesBound / 2)) {
        undist = (totalSamplesBound - 2 - undist);
      }
      undistinguishedCounts[i][distinguished] = undist;
      // cout << "result: " << undist << endl;
      // cout << i << "\t" << distinguished << "\t" << undist << endl;
    }
  }
  cout << " Done.\n" << endl;
}

unsigned int DataFastSMC::readMap(string inFileRoot)
{
  FileUtils::AutoGzIfstream hapsBr;
  if (FileUtils::fileExists(inFileRoot + ".map.gz")) {
    hapsBr.openOrExit(inFileRoot + ".map.gz");
  } else if (FileUtils::fileExists(inFileRoot + ".map")) {
    hapsBr.openOrExit(inFileRoot + ".map");
  } else {
    cerr << "ERROR. Could not find hap file in " + inFileRoot + ".map.gz or " + inFileRoot + ".map" << endl;
    exit(1);
  }
  string line;
  unsigned int pos = 0;
  geneticPositions = vector<double>(sites);
  recRateAtMarker = vector<float>(sites);
  physicalPositions = vector<unsigned long int>(sites);
  while (getline(hapsBr, line)) {
    vector<string> splitStr;
    istringstream iss(line);
    string buf;
    while (iss >> buf)
      splitStr.push_back(buf);
    string ID = splitStr[1];
    float gen = StringUtils::stof(splitStr[2]) / 100.f;
    unsigned long int phys = std::stoul(splitStr[3]);
    geneticPositions[pos] = gen;
    physicalPositions[pos] = phys;
    physicalPositionsMap[phys] = pos;
    if (pos > 0) {
      float genDistFromPrevious = geneticPositions[pos] - geneticPositions[pos - 1];
      unsigned long int physDistFromPrevious = physicalPositions[pos] - physicalPositions[pos - 1];
      float recRate = genDistFromPrevious / physDistFromPrevious;
      recRateAtMarker[pos] = recRate;
      if (pos == 1) {
        // if it's first, add it again for marker 0. Using rate to next marker instead of previous marker
        recRateAtMarker[pos - 1] = recRate;
      }
    }
    pos++;
  }
  cout << "Read " << pos << " markers from map file." << endl;
  if (pos != sites) {
    cerr << "ERROR. Read " << pos << " from map file, expected " << sites << endl;
    exit(1);
  }
  return pos;
}

void DataFastSMC::readSamplesList(string inFileRoot, int jobID, int jobs)
{
  string line;
  // Read samples file
  FileUtils::AutoGzIfstream bufferedReader;
  if (FileUtils::fileExists(inFileRoot + ".samples")) {
    bufferedReader.openOrExit(inFileRoot + ".samples");
  } else if (FileUtils::fileExists(inFileRoot + ".sample")) {
    bufferedReader.openOrExit(inFileRoot + ".sample");
  } else {
    cerr << "ERROR. Could not find sample file in " + inFileRoot + ".sample or " + inFileRoot + ".samples" << endl;
    exit(1);
  }

  unsigned int cpt = 0; // number of current line being processed
  // BufferedReader bufferedReader = Utils.openFile(inFileRoot + ".samples");
  while (getline(bufferedReader, line)) {
    vector<string> splitStr;
    istringstream iss(line);
    string buf;
    while (iss >> buf)
      splitStr.push_back(buf);
    // Skip first two lines (header) if present
    if ((splitStr[0] == "ID_1" && splitStr[1] == "ID_2" && splitStr[2] == "missing") ||
        (splitStr[0] == "0" && splitStr[1] == "0" && splitStr[2] == "0")) {
      continue;
    }
    if ((cpt >= (uint)((w_i - 1) * windowSize) / 2 && cpt < (uint)(w_i * windowSize) / 2) ||
        (cpt >= (uint)((w_j - 1) * windowSize) / 2 && cpt < (uint)(w_j * windowSize) / 2) ||
        (jobs == jobID && cpt >= (uint)((w_j - 1) * windowSize) / 2)) {
      string famId = splitStr[0];
      string IId = splitStr[1];
      FamIDList.push_back(famId);
      IIDList.push_back(IId);
    }
    cpt++;
  }
  bufferedReader.close();
  cout << "Read " << cpt << " samples." << endl;
}

unsigned int DataFastSMC::countSamplesLines(string inFileRoot)
{
  FileUtils::AutoGzIfstream sampleBr;
  if (FileUtils::fileExists(inFileRoot + ".samples")) {
    sampleBr.openOrExit(inFileRoot + ".samples");
  } else if (FileUtils::fileExists(inFileRoot + ".sample")) {
    sampleBr.openOrExit(inFileRoot + ".sample");
  } else {
    cerr << "ERROR. Could not find sample file in " + inFileRoot + ".sample or " + inFileRoot + ".samples" << endl;
    exit(1);
  }
  string line;
  unsigned int cpt = 0;
  while (getline(sampleBr, line)) {
    vector<string> splitStr;
    istringstream iss(line);
    string buf;
    while (iss >> buf)
      splitStr.push_back(buf);
    // Skip first two lines (header) if present
    if ((splitStr[0] == "ID_1" && splitStr[1] == "ID_2" && splitStr[2] == "missing") ||
        (splitStr[0] == "0" && splitStr[1] == "0" && splitStr[2] == "0")) {
      continue;
    }

    cpt++;
  }
  sampleBr.close();
  return cpt;
}

unsigned int DataFastSMC::countHapLines(string inFileRoot)
{
  FileUtils::AutoGzIfstream hapsBr;
  if (FileUtils::fileExists(inFileRoot + ".hap.gz")) {
    hapsBr.openOrExit(inFileRoot + ".hap.gz");
  } else if (FileUtils::fileExists(inFileRoot + ".hap")) {
    hapsBr.openOrExit(inFileRoot + ".hap");
  } else if (FileUtils::fileExists(inFileRoot + ".haps.gz")) {
    hapsBr.openOrExit(inFileRoot + ".haps.gz");
  } else if (FileUtils::fileExists(inFileRoot + ".haps")) {
    hapsBr.openOrExit(inFileRoot + ".haps");
  } else {
    cerr << "ERROR. Could not find hap file in " + inFileRoot + ".hap.gz, " + inFileRoot + ".hap, " + ".haps.gz, or " +
                inFileRoot + ".haps"
         << endl;
    exit(1);
  }
  string line;
  unsigned int pos = 0;
  while (getline(hapsBr, line)) {
    pos++;
  }
  return pos;
}

void DataFastSMC::readHaps(string inFileRoot, bool foldToMinorAlleles, int jobID, int jobs,
                           vector<pair<unsigned long int, double>>& genetic_map)
{
  FileUtils::AutoGzIfstream hapsBr;
  if (FileUtils::fileExists(inFileRoot + ".hap.gz")) {
    hapsBr.openOrExit(inFileRoot + ".hap.gz");
  } else if (FileUtils::fileExists(inFileRoot + ".hap")) {
    hapsBr.openOrExit(inFileRoot + ".hap");
  } else if (FileUtils::fileExists(inFileRoot + ".haps.gz")) {
    hapsBr.openOrExit(inFileRoot + ".haps.gz");
  } else if (FileUtils::fileExists(inFileRoot + ".haps")) {
    hapsBr.openOrExit(inFileRoot + ".haps");
  } else {
    cerr << "ERROR. Could not find hap file in " + inFileRoot + ".hap.gz, " + inFileRoot + ".hap, " + ".haps.gz, or " +
                inFileRoot + ".haps"
         << endl;
    exit(1);
  }
  string line;
  unsigned int pos = 0, monomorphic = 0;
  totalSamplesCount = vector<unsigned int>(sites);
  derivedAlleleCounts = vector<unsigned int>(sites);
  string chr, snpID;
  unsigned long int bp;
  string alleleA, alleleB;
  unsigned int cur_g = 0;
  while (hapsBr >> chr >> snpID >> bp >> alleleA >> alleleB) {
    getline(hapsBr, line);
    int DAcount = 0;
    if (!(line.length() == 4 * sampleSize || line.length() == 4 * sampleSize + 1)) {
      cerr << "ERROR: haps line has wrong length." << endl;
      exit(1);
    }
    if (cur_g == 0) {
      // splitting inputs
      std::string delimiter = ":";
      size_t chrpos = chr.find(delimiter);
      if (std::string::npos == chrpos) {
        chrNumber = std::stoi(chr);
      } else {
        chrNumber = std::stoi(chr.substr(0, chrpos));
      }
      // check if chrNumber is positive and smaller than 1260, which is the max nb of chromosomes (Ophioglossum)
      // if it is not, set a default value to 0
      if (chrNumber <= 0 || chrNumber > 1260) {
        chrNumber = 0;
      }
    }
    readGeneticMap(bp, genetic_map, cur_g, pos);
    int totalSamples = 2 * sampleSize;
    totalSamplesCount[pos] = totalSamples;
    for (uint i = 0; i < 2 * sampleSize; i++) {
      if (line[2 * i + 1] == '1') {
        DAcount++;
      }
    }
    bool minorAlleleValue = foldToMinorAlleles ? (DAcount <= totalSamples - DAcount) : true;
    siteWasFlippedDuringFolding[pos] = !minorAlleleValue;
    // cout << (DAcount <= totalSamples - DAcount) << "\t" << siteWasFlippedDuringFolding[pos] << endl;
    uint cpt = 0;
    for (uint d = 0; d < sampleSize; d++) {
      if ((d >= (uint)((w_i - 1) * windowSize) / 2 && d < (uint)(w_i * windowSize) / 2) ||
          (d >= (uint)((w_j - 1) * windowSize) / 2 && d < (uint)(w_j * windowSize) / 2) ||
          (jobs == jobID && d >= (uint)((w_j - 1) * windowSize) / 2)) {
        Individual& ind = individuals[cpt];
        cpt++;
        uint hap1 = 2 * d;
        uint hap2 = 2 * d + 1;

        if (line[2 * hap1 + 1] == '1') {
          ind.setGenotype(hap1 % 2 + 1, pos, minorAlleleValue);
        } else if (line[2 * hap1 + 1] == '0') {
          ind.setGenotype(hap1 % 2 + 1, pos, !minorAlleleValue);
        } else {
          cerr << "ERROR: hap is not '0' or '1'" << endl;
          exit(1);
        }

        if (line[2 * hap2 + 1] == '1') {
          ind.setGenotype(hap2 % 2 + 1, pos, minorAlleleValue);
        } else if (line[2 * hap2 + 1] == '0') {
          ind.setGenotype(hap2 % 2 + 1, pos, !minorAlleleValue);
        } else {
          cerr << "ERROR: hap is not '0' or '1'" << endl;
          exit(1);
        }
      }
    }

    if (foldToMinorAlleles) {
      derivedAlleleCounts[pos] = std::min(DAcount, totalSamples - DAcount);
    } else {
      derivedAlleleCounts[pos] = DAcount;
    }

    if (DAcount == 0 || DAcount == totalSamples) {
      monomorphic++;
    }
    pos++;
  }
  hapsBr.close();
  cout << "Read data for " << sampleSize * 2 << " haploid samples and " << pos << " markers, " << monomorphic
       << " of which are monomorphic. This job will focus on " << FamIDList.size() * 2 << " haploid samples." << endl;
}

void DataFastSMC::readGeneticMap(unsigned long int bp, vector<pair<unsigned long int, double>>& genetic_map,
                                 unsigned int& cur_g, unsigned int pos)
{
  double cm;
  while (bp > genetic_map[cur_g].first && cur_g < genetic_map.size() - 1)
    cur_g++;
  if (bp >= genetic_map[cur_g].first) {
    // we found this exact marker, or we reached the end of the map
    cm = genetic_map[cur_g].second;
    addMarker(bp, cm, pos);
  } else if (cur_g == 0) {
    // if we haven't hit the map yet, store first map entry
    cm = genetic_map[cur_g].second;
    addMarker(bp, cm, pos);
  } else {
    // interpolate from previous marker
    cm = genetic_map[cur_g - 1].second + (bp - genetic_map[cur_g - 1].first) *
                                             (genetic_map[cur_g].second - genetic_map[cur_g - 1].second) /
                                             (genetic_map[cur_g].first - genetic_map[cur_g - 1].first);
    addMarker(bp, cm, pos);
  }
}

void DataFastSMC::addMarker(unsigned long int physicalPos, double geneticPos, unsigned int pos)
{

  geneticPositions.push_back(geneticPos / 100.f);
  physicalPositions.push_back(physicalPos);
  physicalPositionsMap[physicalPos] = pos;
  if (pos > 0) {
    double genDistFromPrevious = geneticPositions[pos] - geneticPositions[pos - 1];
    unsigned long int physDistFromPrevious = physicalPositions[pos] - physicalPositions[pos - 1];
    float recRate = genDistFromPrevious / physDistFromPrevious;
    if (pos == 1) {
      // if it's first, add it again for marker 0. Using rate to next marker instead of previous marker
      recRateAtMarker.push_back(recRate);
    }
    recRateAtMarker.push_back(recRate);
  }
}
