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
  std::vector<std::string> FamIDList = {};
  std::vector<std::string> IIDList = {};
  std::vector<std::string> famAndIndNameList = {};
  std::vector<Individual> individuals = {};

  unsigned long sampleSize = 0ul;
  unsigned long haploidSampleSize = 0ul;
  int sites = 0;
  int totalSamplesBound = 0;
  bool decodingUsesCSFS = false;
  bool mJobbing = false;
  std::vector<float> geneticPositions = {};
  std::vector<int> physicalPositions = {};
  std::vector<bool> siteWasFlippedDuringFolding = {};
  std::vector<float> recRateAtMarker = {};
  std::vector<std::vector<int>> undistinguishedCounts = {};


  // Variables relating to FastSMC
  int chrNumber = 0;
  unsigned int windowSize = 0u; // window size in triangles for each job
  unsigned int w_i = 0u;        // window id for ind_i for jobs
  unsigned int w_j = 0u;        // window id for ind_j for jobs
  bool is_j_above_diag = false;
  std::unordered_map<int, unsigned int> physicalPositionsMap = {}; // map where key=physicalPosition, value=indexPosition

  /**
   * Construct the data object
   *
   * @param inFileRoot location of input files
   * @param totalSamplesBound
   * @param foldToMinorAlleles
   * @param decodingUsesCSFS whether to decode using CSFS
   * @param jobID the jobID which defaults to -1 indicating no jobbing
   * @param jobs the number of jobs which defaults to -1 indicating no jobbing
   * @param useKnownSeed use a known random seed ensuring predictable random numbers for testing
   */
  Data(const std::string& inFileRoot, int totalSamplesBound, bool foldToMinorAlleles,
       bool decodingUsesCSFS, int jobID = -1, int jobs = -1, bool useKnownSeed = false);

  static int countHapLines(std::string inFileRoot);
  static int countSamplesLines(std::string inFileRoot);

private:

  /**
   * TODO: harmonise functionality in this class so that this switch is not necessary
   */
  bool TODO_REMOVE_FASTSMC = false;

  /**
   * Determine whether a sample should be read, based on the jobID, number of jobs, and the number of lines processed.
   * ASMC will always return true, but FastSMC will determine whether to read a sample.
   *
   * @param linesProcessed the number of lines processed so far
   * @param jobID the jobID, which will be the default value of -1 for ASMC
   * @param jobs the number of jobs, which will be the default value of -1 for ASMC
   * @return whether to read the sample
   */
  bool readSample(unsigned linesProcessed, int jobID, int jobs);

  /**
   * Read the samples file and populate members `FamIDList`, `IIDList` and `famAndIndNameList`.
   *
   * @param inFileRoot location of input files
   * @param jobID the jobID which defaults to -1 indicating no jobbing
   * @param jobs the number of jobs which defaults to -1 indicating no jobbing
   */
  void readSamplesList(const std::string& inFileRoot, int jobID, int jobs);

  void readHaps(std::string inFileRoot, bool foldToMinorAlleles);
  void readHaps(std::string inFileRoot, bool foldToMinorAlleles, int jobID, int jobs,
                std::vector<std::pair<unsigned long int, double>>& genetic_map);

  int readMap(std::string inFileRoot);

  /**
   * Subsumed functionality from FastSMC to read genetic map as a vector of pairs.
   * TODO: can this be harmonised with the other readMap method?
   * @param inFileRoot
   * @return
   */
  static std::vector<std::pair<unsigned long, double>> readMapFastSMC(const std::string& inFileRoot);


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
