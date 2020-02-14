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

#ifndef ASMC_DATA_HPP
#define ASMC_DATA_HPP

#include <string>
#include <vector>

#include "Individual.hpp"

class Data
{

public:
  std::vector<std::string> FamIDList;
  std::vector<std::string> IIDList;
  std::vector<std::string> famAndIndNameList;
  std::vector<Individual> individuals;

  int sampleSize;
  int haploidSampleSize;
  int sites;
  int totalSamplesBound;
  bool decodingUsesCSFS = false;
  std::vector<float> geneticPositions;
  std::vector<int> physicalPositions;
  std::vector<bool> siteWasFlippedDuringFolding;
  std::vector<float> recRateAtMarker;
  std::vector<std::vector<int>> undistinguishedCounts;

  Data(std::string hapsFileRoot, int numOfSites, int totalSamplesBound, bool foldToMinorAlleles, bool decodingUsesCSFS);
  static int countHapLines(std::string hapsFileRoot);

private:
  void readSamplesList(std::string hapsFileRoot);
  void readHaps(std::string hapsFileRoot, bool foldToMinorAlleles);
  int readMap(std::string hapsFileRoot);
  void makeUndistinguished(bool foldToMinorAlleles);
  std::vector<int> totalSamplesCount;
  std::vector<int> derivedAlleleCounts;
  std::vector<std::string> SNP_IDs;
  int sampleHypergeometric(int populationSize, int numberOfSuccesses, int sampleSize);
};

#endif // ASMC_DATA_HPP
