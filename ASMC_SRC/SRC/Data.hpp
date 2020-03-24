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
#include <unordered_map>
#include <utility>
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


  // Variables relating to FastSMC
  int chrNumber;
  unsigned int windowSize; // window size in triangles for each job
  unsigned int w_i;        // window id for ind_i for jobs
  unsigned int w_j;        // window id for ind_j for jobs
  bool is_j_above_diag;
  std::unordered_map<int, unsigned int> physicalPositionsMap; // map where key=physicalPosition, value=indexPosition


  Data(std::string inFileRoot, int numOfSites, int totalSamplesBound, bool foldToMinorAlleles, bool decodingUsesCSFS);
  Data(std::string inFileRoot, unsigned int numOfSites, int totalSamplesBound, bool foldToMinorAlleles,
       bool decodingUsesCSFS, int jobID, int jobs, std::vector<std::pair<unsigned long int, double>>& genetic_map);

  static int countHapLines(std::string inFileRoot);
  static int countSamplesLines(std::string inFileRoot);



private:

  void readSamplesList(std::string inFileRoot);
  void readSamplesList(std::string inFileRoot, int jobID, int jobs);


  void readHaps(std::string inFileRoot, bool foldToMinorAlleles);
  void readHaps(std::string inFileRoot, bool foldToMinorAlleles, int jobID, int jobs,
                std::vector<std::pair<unsigned long int, double>>& genetic_map);

  int readMap(std::string inFileRoot);
  void makeUndistinguished(bool foldToMinorAlleles);
  std::vector<int> totalSamplesCount;
  std::vector<int> derivedAlleleCounts;
  std::vector<std::string> SNP_IDs;
  int sampleHypergeometric(int populationSize, int numberOfSuccesses, int sampleSize);


  void readGeneticMap(unsigned long int bp, std::vector<std::pair<unsigned long int, double>>& genetic_map,
                      unsigned int& cur_g, unsigned int pos);

  void addMarker(unsigned long int physicalPosition, double geneticPosition, unsigned int pos);




};

#endif // ASMC_DATA_HPP
