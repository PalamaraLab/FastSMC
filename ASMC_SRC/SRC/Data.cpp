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

inline bool fileExists(const std::string& name)
{
  ifstream f(name.c_str());
  return f.good();
}

Data::Data(string hapsFileRoot, int _sites, int _totalSamplesBound,
    bool foldToMinorAlleles, bool _decodingUsesCSFS)
    : sites(_sites)
    , totalSamplesBound(_totalSamplesBound)
    , decodingUsesCSFS(_decodingUsesCSFS)
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

// unoptimized sampling of hypergeometric
int Data::sampleHypergeometric(
    int populationSize, int numberOfSuccesses, int sampleSize)
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
           << " non-missing individuals, but the CSFS requires " << totalSamplesBound
           << endl;
      exit(1);
    }
    int ancestralAlleles = totalSamples - derivedAlleles;
    if (foldToMinorAlleles && derivedAlleles > ancestralAlleles) {
      cerr << "Minor alleles has frequency > 50%. Data is supposed to be folded.\n";
    }
    for (int distinguished = 0; distinguished < 3; distinguished++) {
      // hypergeometric with (derivedAlleles - distinguished) derived alleles, (samples
      // - 2) samples
      int undist = sampleHypergeometric(
          totalSamples - 2, derivedAlleles - distinguished, totalSamplesBound - 2);
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
  if (fileExists(hapsFileRoot + ".map.gz")) {
    hapsBr.openOrExit(hapsFileRoot + ".map.gz");
  } else if (fileExists(hapsFileRoot + ".map")) {
    hapsBr.openOrExit(hapsFileRoot + ".map");
  } else {
    cerr << "ERROR. Could not find hap file in " + hapsFileRoot + ".map.gz or "
            + hapsFileRoot + ".map"
         << endl;
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
  if (fileExists(hapsFileRoot + ".samples")) {
    bufferedReader.openOrExit(hapsFileRoot + ".samples");
  } else if (fileExists(hapsFileRoot + ".sample")) {
    bufferedReader.openOrExit(hapsFileRoot + ".sample");
  } else {
    cerr << "ERROR. Could not find sample file in " + hapsFileRoot + ".sample or "
            + hapsFileRoot + ".samples"
         << endl;
    exit(1);
  }

  while (getline(bufferedReader, line)) {
    vector<string> splitStr;
    istringstream iss(line);
    string buf;
    while (iss >> buf)
      splitStr.push_back(buf);
    // Skip first two lines (header) if present
    if ((splitStr[0] == "ID_1" && splitStr[1] == "ID_2" && splitStr[2] == "missing")
        || (splitStr[0] == "0" && splitStr[1] == "0" && splitStr[2] == "0")) {
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

int Data::countHapLines(string hapsFileRoot)
{
  FileUtils::AutoGzIfstream hapsBr;
  if (fileExists(hapsFileRoot + ".hap.gz")) {
    hapsBr.openOrExit(hapsFileRoot + ".hap.gz");
  } else if (fileExists(hapsFileRoot + ".hap")) {
    hapsBr.openOrExit(hapsFileRoot + ".hap");
  } else if (fileExists(hapsFileRoot + ".haps.gz")) {
    hapsBr.openOrExit(hapsFileRoot + ".haps.gz");
  } else if (fileExists(hapsFileRoot + ".haps")) {
    hapsBr.openOrExit(hapsFileRoot + ".haps");
  } else {
    cerr << "ERROR. Could not find hap file in " + hapsFileRoot + ".hap.gz, "
            + hapsFileRoot + ".hap, " + ".haps.gz, or " + hapsFileRoot + ".haps"
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

void Data::readHaps(string hapsFileRoot, bool foldToMinorAlleles)
{
  FileUtils::AutoGzIfstream hapsBr;
  if (fileExists(hapsFileRoot + ".hap.gz")) {
    hapsBr.openOrExit(hapsFileRoot + ".hap.gz");
  } else if (fileExists(hapsFileRoot + ".hap")) {
    hapsBr.openOrExit(hapsFileRoot + ".hap");
  } else if (fileExists(hapsFileRoot + ".haps.gz")) {
    hapsBr.openOrExit(hapsFileRoot + ".haps.gz");
  } else if (fileExists(hapsFileRoot + ".haps")) {
    hapsBr.openOrExit(hapsFileRoot + ".haps");
  } else {
    cerr << "ERROR. Could not find hap file in " + hapsFileRoot + ".hap.gz, "
            + hapsFileRoot + ".hap, " + ".haps.gz, or " + hapsFileRoot + ".haps"
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
    if (!(line.length() == 4 * famAndIndNameList.size()
            || line.length() == 4 * famAndIndNameList.size() + 1)) {
      cerr << "ERROR: haps line has wrong length. Length is " << line.length()
           << ", should be 4*" << famAndIndNameList.size() << " = "
           << 4 * famAndIndNameList.size() << "." << endl;
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
    bool minorAlleleValue
        = foldToMinorAlleles ? (DAcount <= totalSamples - DAcount) : true;
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
  cout << "Read data for " << famAndIndNameList.size() * 2 << " haploid samples and "
       << pos << " markers, " << monomorphic << " of which are monomorphic." << endl;
}
