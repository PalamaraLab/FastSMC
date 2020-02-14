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

#ifndef ASMC_DATA_FASTSMC_HPP
#define ASMC_DATA_FASTSMC_HPP

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Individual.hpp"

class DataFastSMC
{

  /* **************************** */
  /* **************************** */
  // loads haps/samples
  /* **************************** */
  /* **************************** */

public:
  /* **************************** */
  // structures for sample names and individuals
  std::vector<std::string> FamIDList;
  std::vector<std::string> IIDList;
  std::vector<Individual> individuals;

  // some variables about the data
  int chrNumber;
  unsigned int sampleSize; // nb of diploid ind
  unsigned int sites;      // nb of sites
  unsigned int totalSamplesBound;
  unsigned int windowSize; // window size in triangles for each job
  unsigned int w_i;        // window id for ind_i for jobs
  unsigned int w_j;        // window id for ind_j for jobs
  bool is_j_above_diag;
  bool decodingUsesCSFS = false;
  std::vector<double>
      geneticPositions; // vector of genetic positions for each site (3rd column of map file divided by 100)
  std::vector<unsigned long int> physicalPositions; // vector of pysical positions for each site
  std::vector<bool> siteWasFlippedDuringFolding;
  std::vector<float> recRateAtMarker; // rec rate for each site
  std::vector<std::vector<int>> undistinguishedCounts;
  std::unordered_map<long int, unsigned int> physicalPositionsMap; // map where key=physicalPosition,
                                                                   // value=indexPosition

  DataFastSMC(std::string inFileRoot, unsigned int numOfSites, int totalSamplesBound, bool foldToMinorAlleles,
              bool decodingUsesCSFS, int jobID, int jobs,
              std::vector<std::pair<unsigned long int, double>>& genetic_map);
  static unsigned int countHapLines(std::string inFileRoot);
  static unsigned int countSamplesLines(std::string inFileRoot);

private:
  std::vector<int> prefixArray(int k);
  void readSamplesList(std::string inFileRoot, int jobID, int jobs);
  void readHaps(std::string inFileRoot, bool foldToMinorAlleles, int jobID, int jobs,
                std::vector<std::pair<unsigned long int, double>>& genetic_map);
  unsigned int readMap(std::string inFileRoot);
  void makeUndistinguished(bool foldToMinorAlleles);
  std::vector<unsigned int> totalSamplesCount;   // vector of nbhaplotypes for each site
  std::vector<unsigned int> derivedAlleleCounts; // vector ... for each site
  int sampleHypergeometric(int populationSize, int numberOfSuccesses, unsigned int sampleSize);
  void readGeneticMap(unsigned long int bp, std::vector<std::pair<unsigned long int, double>>& genetic_map,
                      unsigned int& cur_g, unsigned int pos);
  void addMarker(unsigned long int physicalPosition, double geneticPosition, unsigned int pos);
};

#endif // ASMC_DATA_FASTSMC_HPP
