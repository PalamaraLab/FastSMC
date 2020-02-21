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

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "FileUtils.hpp"
#include "StringUtils.hpp"
#include "Types.hpp"

#include "Data.hpp"
#include <boost/math/distributions/hypergeometric.hpp>

using namespace std;

Data::Data(string hapsFileRoot, int _sites, int _totalSamplesBound, bool foldToMinorAlleles, bool _decodingUsesCSFS)
    : sites(_sites), totalSamplesBound(_totalSamplesBound), decodingUsesCSFS(_decodingUsesCSFS)
{
  readSamplesList(hapsFileRoot);
  // now we have the sample sizes
  sampleSize = static_cast<int>(famAndIndNameList.size());
  haploidSampleSize = sampleSize * 2;
  siteWasFlippedDuringFolding = vector<bool>(sites, false);
  // and the number of sites
  for (uint i = 0; i < famAndIndNameList.size(); i++) {
    individuals.push_back(Individual(sites));
  }
  // now read all the data
  readHaps(hapsFileRoot, foldToMinorAlleles);
  readMap(hapsFileRoot);
  // make undistinguished counts
  if (decodingUsesCSFS) {
    makeUndistinguished(foldToMinorAlleles);
  }
}

Data::Data(string inFileRoot, unsigned int _sites, int _totalSamplesBound, bool foldToMinorAlleles,
           bool _decodingUsesCSFS, int jobID, int jobs, vector<pair<unsigned long int, double>>& genetic_map)
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

// unoptimized sampling of hypergeometric
int Data::sampleHypergeometric(int populationSize, int numberOfSuccesses, int sampleSize)
{
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
  for (int i = 0; i < sampleSize; i++) {
    ret += samplingVector[i];
  }
  return ret;
}

void Data::makeUndistinguished(bool foldToMinorAlleles)
{
  // cout << "Building undistinguished counts...";
  undistinguishedCounts = vector<vector<int>>(derivedAlleleCounts.size());
  for (uint i = 0; i < derivedAlleleCounts.size(); i++) {
    undistinguishedCounts[i] = vector<int>(3);
    int derivedAlleles = derivedAlleleCounts[i];
    int totalSamples = totalSamplesCount[i];
    if (decodingUsesCSFS && totalSamplesBound > totalSamples) {
      cerr << "ERROR. SNP " << SNP_IDs[i] << " has " << totalSamples
           << " non-missing individuals, but the CSFS requires " << totalSamplesBound << endl;
      exit(1);
    }
    int ancestralAlleles = totalSamples - derivedAlleles;
    if (foldToMinorAlleles && derivedAlleles > ancestralAlleles) {
      cerr << "Minor alleles has frequency > 50%. Data is supposed to be folded.\n";
    }
    for (int distinguished = 0; distinguished < 3; distinguished++) {
      // hypergeometric with (derivedAlleles - distinguished) derived alleles, (samples
      // - 2) samples
      int undist = sampleHypergeometric(totalSamples - 2, derivedAlleles - distinguished, totalSamplesBound - 2);
      if (foldToMinorAlleles && (undist + distinguished > totalSamplesBound / 2)) {
        undist = (totalSamplesBound - 2 - undist);
      }
      undistinguishedCounts[i][distinguished] = undist;
    }
  }
  // cout << " Done." << endl;
}

int Data::readMap(string hapsFileRoot)
{
  FileUtils::AutoGzIfstream hapsBr;
  if (FileUtils::fileExists(hapsFileRoot + ".map.gz")) {
    hapsBr.openOrExit(hapsFileRoot + ".map.gz");
  } else if (FileUtils::fileExists(hapsFileRoot + ".map")) {
    hapsBr.openOrExit(hapsFileRoot + ".map");
  } else {
    cerr << "ERROR. Could not find hap file in " + hapsFileRoot + ".map.gz or " + hapsFileRoot + ".map" << endl;
    exit(1);
  }
  string line;
  int pos = 0;
  SNP_IDs = vector<string>(sites);
  geneticPositions = vector<float>(sites);
  recRateAtMarker = vector<float>(sites);
  physicalPositions = vector<int>(sites);
  while (getline(hapsBr, line)) {
    vector<string> splitStr;
    istringstream iss(line);
    string buf;
    while (iss >> buf)
      splitStr.push_back(buf);
    string ID = splitStr[1];
    float gen = StringUtils::stof(splitStr[2]) / 100.f;
    int phys = std::stoi(splitStr[3]);
    SNP_IDs[pos] = ID;
    geneticPositions[pos] = gen;
    physicalPositions[pos] = phys;
    if (pos > 0) {
      float genDistFromPrevious = geneticPositions[pos] - geneticPositions[pos - 1];
      int physDistFromPrevious = physicalPositions[pos] - physicalPositions[pos - 1];
      float recRate = genDistFromPrevious / physDistFromPrevious;
      recRateAtMarker[pos] = recRate;
      if (pos == 1) {
        // if it's first, add it again for marker 0. Using rate to next marker instead
        // of previous marker
        recRateAtMarker[pos] = recRate;
      }
    }
    pos++;
  }
  // cout << "\tRead " << pos << " markers from map file." << endl;
  if (pos != sites) {
    cerr << "ERROR. Read " << pos << " from map file, expected " << sites << endl;
    exit(1);
  }
  return pos;
}

void Data::readSamplesList(string hapsFileRoot)
{
  string line;
  // Read samples file
  FileUtils::AutoGzIfstream bufferedReader;
  if (FileUtils::fileExists(hapsFileRoot + ".samples")) {
    bufferedReader.openOrExit(hapsFileRoot + ".samples");
  } else if (FileUtils::fileExists(hapsFileRoot + ".sample")) {
    bufferedReader.openOrExit(hapsFileRoot + ".sample");
  } else {
    cerr << "ERROR. Could not find sample file in " + hapsFileRoot + ".sample or " + hapsFileRoot + ".samples" << endl;
    exit(1);
  }

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
    string famId = splitStr[0];
    string IId = splitStr[1];
    string name = famId + "\t" + IId;
    FamIDList.push_back(famId);
    IIDList.push_back(IId);
    famAndIndNameList.push_back(name);
  }
  bufferedReader.close();
  // cout << "\tRead " << famAndIndNameList.size() << " samples." << endl;
}

void Data::readSamplesList(string inFileRoot, int jobID, int jobs)
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

int Data::countHapLines(string hapsFileRoot)
{
  FileUtils::AutoGzIfstream hapsBr;
  if (FileUtils::fileExists(hapsFileRoot + ".hap.gz")) {
    hapsBr.openOrExit(hapsFileRoot + ".hap.gz");
  } else if (FileUtils::fileExists(hapsFileRoot + ".hap")) {
    hapsBr.openOrExit(hapsFileRoot + ".hap");
  } else if (FileUtils::fileExists(hapsFileRoot + ".haps.gz")) {
    hapsBr.openOrExit(hapsFileRoot + ".haps.gz");
  } else if (FileUtils::fileExists(hapsFileRoot + ".haps")) {
    hapsBr.openOrExit(hapsFileRoot + ".haps");
  } else {
    cerr << "ERROR. Could not find hap file in " + hapsFileRoot + ".hap.gz, " + hapsFileRoot + ".hap, " +
                ".haps.gz, or " + hapsFileRoot + ".haps"
         << endl;
    exit(1);
  }
  string line;
  int pos = 0;
  while (getline(hapsBr, line)) {
    pos++;
  }
  return pos;
}

int Data::countSamplesLines(string inFileRoot)
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
  int cpt = 0;
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

void Data::readHaps(string hapsFileRoot, bool foldToMinorAlleles)
{
  FileUtils::AutoGzIfstream hapsBr;
  if (FileUtils::fileExists(hapsFileRoot + ".hap.gz")) {
    hapsBr.openOrExit(hapsFileRoot + ".hap.gz");
  } else if (FileUtils::fileExists(hapsFileRoot + ".hap")) {
    hapsBr.openOrExit(hapsFileRoot + ".hap");
  } else if (FileUtils::fileExists(hapsFileRoot + ".haps.gz")) {
    hapsBr.openOrExit(hapsFileRoot + ".haps.gz");
  } else if (FileUtils::fileExists(hapsFileRoot + ".haps")) {
    hapsBr.openOrExit(hapsFileRoot + ".haps");
  } else {
    cerr << "ERROR. Could not find hap file in " + hapsFileRoot + ".hap.gz, " + hapsFileRoot + ".hap, " +
                ".haps.gz, or " + hapsFileRoot + ".haps"
         << endl;
    exit(1);
  }
  string line;
  int pos = 0, monomorphic = 0;
  totalSamplesCount = vector<int>(sites);
  derivedAlleleCounts = vector<int>(sites);
  string chr, snpID;
  int bp;
  string alleleA, alleleB;
  while (hapsBr >> chr >> snpID >> bp >> alleleA >> alleleB) {
    getline(hapsBr, line);
    int DAcount = 0;
    if (!(line.length() == 4 * famAndIndNameList.size() || line.length() == 4 * famAndIndNameList.size() + 1)) {
      cerr << "ERROR: haps line has wrong length. Length is " << line.length() << ", should be 4*"
           << famAndIndNameList.size() << " = " << 4 * famAndIndNameList.size() << "." << endl;
      cerr << "\thaps line is: " << line << endl;

      exit(1);
    }
    int totalSamples = static_cast<int>(2 * famAndIndNameList.size());
    totalSamplesCount[pos] = totalSamples;
    for (uint i = 0; i < 2 * famAndIndNameList.size(); i++) {
      if (line[2 * i + 1] == '1') {
        DAcount++;
      }
    }
    bool minorAlleleValue = foldToMinorAlleles ? (DAcount <= totalSamples - DAcount) : true;
    siteWasFlippedDuringFolding[pos] = !minorAlleleValue;
    for (uint i = 0; i < 2 * famAndIndNameList.size(); i++) {
      int indIndex = i / 2;
      Individual& ind = individuals[indIndex];
      if (line[2 * i + 1] == '1') {
        ind.setGenotype(i % 2 + 1, pos, minorAlleleValue);
      } else if (line[2 * i + 1] == '0') {
        ind.setGenotype(i % 2 + 1, pos, !minorAlleleValue);
      } else {
        cerr << "ERROR: hap is not '0' or '1'" << endl;
        exit(1);
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
  cout << "Read data for " << famAndIndNameList.size() * 2 << " haploid samples and " << pos << " markers, "
       << monomorphic << " of which are monomorphic." << endl;
}

void Data::readHaps(string inFileRoot, bool foldToMinorAlleles, int jobID, int jobs,
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
  totalSamplesCount = vector<int>(sites);
  derivedAlleleCounts = vector<int>(sites);
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

void Data::readGeneticMap(unsigned long int bp, vector<pair<unsigned long int, double>>& genetic_map,
                          unsigned int& cur_g, unsigned int pos)
{
  double cm;

  while (bp > genetic_map[cur_g].first && cur_g < genetic_map.size() - 1) {
    cur_g++;
  }

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

void Data::addMarker(unsigned long int physicalPos, double geneticPos, unsigned int pos)
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