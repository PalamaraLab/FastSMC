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


#ifndef DATA_HPP
#define DATA_HPP

#include <vector>
#include <string>

#include "Individual.hpp"

using namespace std;

class Data {

public:

  vector <string> FamIDList;
  vector <string> IIDList;
  vector <string> famAndIndNameList;
  vector <Individual> individuals;

  int sampleSize;
  int haploidSampleSize;
  int sites;
  int totalSamplesBound;
  bool decodingUsesCSFS = false;
  vector<float> geneticPositions;
  vector<int> physicalPositions;
  vector<bool> siteWasFlippedDuringFolding;
  vector<float> recRateAtMarker;
  vector< vector<int> > undistinguishedCounts;

  Data(string hapsFileRoot, int numOfSites, int totalSamplesBound, bool foldToMinorAlleles, bool decodingUsesCSFS);
  static int countHapLines(string hapsFileRoot);

private:
  void readSamplesList(string hapsFileRoot);
  void readHaps(string hapsFileRoot, bool foldToMinorAlleles);
  int readMap(string hapsFileRoot);
  void makeUndistinguished(bool foldToMinorAlleles);
  vector<int> totalSamplesCount;
  vector<int> derivedAlleleCounts;
  vector<string> SNP_IDs;
  int sampleHypergeometric(int populationSize, int numberOfSuccesses, int sampleSize);

};

#endif
