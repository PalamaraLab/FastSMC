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

#include "catch.hpp"

#include "HMM.hpp"
#include "FileUtils.hpp"
#include <sstream>



#include <bitset>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <inttypes.h>
#include <iostream>
#include <list>
#include <math.h>
#include <sstream>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <vector>

// Includes for seeding
#include <boost/unordered_map.hpp>

#include "Data.hpp"
#include "DecodingParams.hpp"
#include "DecodingQuantities.hpp"
#include "FileUtils.hpp"
#include "HMM.hpp"
#include "StringUtils.hpp"
#include "Timer.hpp"

#include "HASHING/ExtendHash.hpp"
#include "HASHING/Individuals.hpp"
#include "HASHING/SeedHash.hpp"

using namespace std;

vector<int> hist_ctr;
vector<float> hist_win;

// Set the word/hash sizes here at compile time, must be able to cast from a ulong
float PAR_MIN_MATCH = 1;
bool PAR_DIAGNOSTICS = false;
int PAR_GAP = 1;
int MAX_seeds = 0;
int GLOBAL_READ_WORDS = 0;
int GLOBAL_CURRENT_WORD = 0;
int GLOBAL_SKIPPED_WORDS = 0;
int GEN_THRESHOLD = 100;
bool PAR_BIN_OUT = false;
bool PAR_HAPLOID = true;

bool is_j_above_diag;
unsigned int windowSize;
unsigned int w_i;
unsigned int w_j;
int jobID, jobs;

vector<float>* all_markers;
unsigned long int num_ind;
unsigned long int num_ind_tot;

const int CONST_READ_AHEAD = 10;
const int WORD_SIZE = 64;
typedef uint64_t hash_size;

struct Marker {
  string id;
  unsigned long int pos;
  double cm;

  string print()
  {
    stringstream ss;
    ss << id << '\t' << pos << '\t' << cm << endl;
    return ss.str();
  }
};

vector<Individuals<WORD_SIZE, CONST_READ_AHEAD>> all_ind;

bool isHapInJob(unsigned int i)
{
  return ((i >= (w_i - 1) * windowSize && i < w_i * windowSize) ||
          (i >= (w_j - 1) * windowSize && i < w_j * windowSize) || (jobs == jobID && i >= (w_j - 1) * windowSize));
}

bool isSampleInJob(unsigned int i)
{
  return ((i >= (uint)((w_i - 1) * windowSize) / 2 && i < (uint)(w_i * windowSize) / 2) ||
          (i >= (uint)((w_j - 1) * windowSize) / 2 && i < (uint)(w_j * windowSize) / 2) ||
          (jobs == jobID && i >= (uint)((w_j - 1) * windowSize) / 2));
}

double get_cpu_time()
{
  return (double)clock() / CLOCKS_PER_SEC;
}




TEST_CASE("test FastSMC HMM with regression test", "[FastSMC_regression]")
{
  // Setup will go here
  srand(1234);

  DecodingParams params;
  params.decodingQuantFile = ASMC_FILE_DIR "/FASTSMC_EXAMPLE/out.25.n300.chr2.len30.dens1.disc10-20-2000.demoCEU.mapnorm.array.decodingQuantities.gz";
  params.inFileRoot = ASMC_FILE_DIR "/FASTSMC_EXAMPLE/out.25.n300.chr2.len30.dens1.disc10-20-2000.demoCEU.mapnorm.array";
  params.outFileRoot = "/tmp/FastSMCresults";
  params.decodingModeString = "array";
  params.decodingMode = DecodingMode::arrayFolded;
  params.foldData = true;
  params.usingCSFS = true;
  params.map = ASMC_FILE_DIR "/FASTSMC_EXAMPLE/chr2.map";
  params.batchSize = 32;
  params.recallThreshold = 3;
  params.min_m = 1.5;
  params.GERMLINE = true;
  params.BIN_OUT = false;
  params.time = 50;
  params.noConditionalAgeEstimates = true;
  params.doPerPairMAP = true;
  params.doPerPairPosteriorMean = true;

  double TIME_start = get_cpu_time();
  double TIME_prev = TIME_start;
  float PAR_MIN_MAF = 0;
  float PAR_skip = 0;

  string line, discard;
  bool opt_error = 0;
  int c;

  PAR_MIN_MATCH = params.min_m;
  PAR_MIN_MAF = params.min_maf;
  PAR_GAP = params.gap;
  PAR_skip = params.skip;
  MAX_seeds = params.max_seeds;

  // read decoding quantities from file
  std::string str(params.decodingQuantFile.c_str());
  DecodingQuantities decodingQuantities(params.decodingQuantFile.c_str());

  unsigned int sequenceLength = Data::countHapLines(params.inFileRoot.c_str());

  FileUtils::AutoGzIfstream file_haps;
  file_haps.openOrExit(params.inFileRoot + ".hap.gz");

  FileUtils::AutoGzIfstream file_samp;
  file_samp.openOrExit(params.inFileRoot + ".samples");

  ifstream file_genm(params.map);

  string map_field[3];
  stringstream ss;

  // *** read genetic map
  vector<pair<unsigned long int, double>> genetic_map;
  unsigned int cur_g = 0;
  while (getline(file_genm, line)) {
    ss.clear();
    ss.str(line);
    ss >> map_field[0] >> map_field[1] >> map_field[2];
    if (map_field[0] == "position" || map_field[0] == "")
      continue;
    genetic_map.push_back(pair<unsigned long int, double>(stol(map_field[0]), stod(map_field[2])));
    cur_g++;
  }
  file_genm.close();

  Data data(params.inFileRoot, sequenceLength, decodingQuantities.CSFSSamples, params.foldData, params.usingCSFS,
            params.jobInd, params.jobs, genetic_map);

  // freeing memory
  genetic_map.clear();
  vector<pair<unsigned long int, double>>().swap(genetic_map);

  HMM hmm(data, decodingQuantities, params, !params.noBatches);
  hmm.decodeAll(params.jobs, params.jobInd);

  if (!params.GERMLINE) {
    exit(1);
  }

  jobs = params.jobs;
  jobID = params.jobInd;
  windowSize = data.windowSize;
  w_i = data.w_i;
  w_j = data.w_j;
  is_j_above_diag = data.is_j_above_diag;

  // *** read Individuals information
  unsigned int linectr = 0;
  while (getline(file_samp, line)) {
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

    ss.clear();
    ss.str(line);
    ss >> map_field[0] >> map_field[1];
    if (isSampleInJob(linectr)) {
      if (PAR_HAPLOID) {
        all_ind.emplace_back(2 * linectr);
        all_ind.emplace_back(2 * linectr + 1);
      } else {
        all_ind.emplace_back(2 * linectr);
        all_ind.emplace_back(2 * linectr);
      }
    }
    linectr++;
  }
  file_samp.close();
  num_ind_tot = data.sampleSize * 2;
  num_ind = all_ind.size();

  Marker cur_marker;
  // track position through genetic map
  cur_g = 0;
  unsigned int snp_ctr;
  char al[2], inp;
  string cur_al;

  // Storage for seeds
  SeedHash<WORD_SIZE, CONST_READ_AHEAD> seeds;

  hash_size word[2];

  // Storage for extensions
  ExtendHash<WORD_SIZE> extend(PAR_HAPLOID ? num_ind : num_ind / 2, PAR_HAPLOID);

  // Hash individual words
  GLOBAL_READ_WORDS = 0;
  GLOBAL_CURRENT_WORD = 0;

  string marker_id;
  unsigned long int marker_pos;
  all_markers = &data.geneticPositions;
  while (1) {
    snp_ctr = 0;
    while (getline(file_haps, line)) {
      // read the meta data
      ss.clear();
      ss.str(line);
      ss >> map_field[0] >> marker_id >> marker_pos >> al[0] >> al[1];
      if (map_field[0] == "")
        continue;

      // restrict on MAF
      if (PAR_MIN_MAF > 0) {
        int maf_ctr = 0;
        for (int i = 0; i < num_ind_tot; i++) {
          ss >> inp;
          if (inp == '1')
            maf_ctr++;
        }
        float maf = (float)maf_ctr / num_ind_tot;
        if (maf < PAR_MIN_MAF || maf > 1 - PAR_MIN_MAF)
          continue;

        // re-load the data
        ss.clear();
        ss.str(line);
        ss >> map_field[0] >> map_field[0] >> map_field[0] >> al[0] >> al[1];
      }

      // read haplotype
      unsigned int hap_ctr = 0;
      for (unsigned int i = 0; i < num_ind_tot; i++) {
        ss >> inp;
        if (isSampleInJob(i / 2)) {
          if (inp == '1')
            all_ind[hap_ctr].setMarker(GLOBAL_READ_WORDS, snp_ctr);
          hap_ctr++;
        }
      }
      snp_ctr++;

      if (snp_ctr % WORD_SIZE == 0) {
        if (++GLOBAL_READ_WORDS >= CONST_READ_AHEAD)
          break;
        else
          cerr << "*** loading word buffer " << GLOBAL_READ_WORDS << " / " << CONST_READ_AHEAD << endl;
        snp_ctr = 0;
      }
    }

    // end if read all data
    if (GLOBAL_CURRENT_WORD >= GLOBAL_READ_WORDS) {
      break;
    }

    for (unsigned int i = 0; i < num_ind; i++) {
      seeds.insertIndividuals(i, all_ind[i].getWordHash(GLOBAL_CURRENT_WORD));
    }

    GLOBAL_SKIPPED_WORDS = 0;
    int cur_seeds = seeds.size();
    unsigned long cur_pairs = 0;

    // skip low-complexity words
    if ((float)cur_seeds / num_ind > PAR_skip) {
      cur_pairs = seeds.extendAllPairs(&extend, GLOBAL_CURRENT_WORD, all_ind, MAX_seeds, jobID, jobs, w_i, w_j,
                                       windowSize, GLOBAL_READ_WORDS, GLOBAL_SKIPPED_WORDS, GLOBAL_CURRENT_WORD,
                                       is_j_above_diag);
      extend.clearPairsPriorTo(GLOBAL_CURRENT_WORD - PAR_GAP, GLOBAL_CURRENT_WORD, PAR_MIN_MATCH, *all_markers, hmm);
    } else {
      cerr << "low complexity word - " << cur_seeds << " - skipping" << endl;
      extend.extendAllPairsTo(GLOBAL_CURRENT_WORD);
    }

    seeds.clear();

    for (unsigned int i = 0; i < num_ind; i++) {
      all_ind[i].clear(GLOBAL_CURRENT_WORD);
    }
    GLOBAL_CURRENT_WORD++;
  }

  extend.clearAllPairs(PAR_MIN_MATCH, *all_markers, hmm);
  file_haps.close();

  hmm.finishFromGERMLINE();

  SECTION("regression test")
  {
    // Read lines from existing regression test output into a vector of strings
    std::vector<std::string> regressionLines;
    regressionLines.reserve(1535);
    {
      FileUtils::AutoGzIfstream fin_regression;
      fin_regression.openOrExit(ASMC_FILE_DIR "/FASTSMC_EXAMPLE/regression_output.ibd.gz");
      for (std::string f_line; getline(fin_regression, f_line);) {
        regressionLines.emplace_back(f_line);
      }
      fin_regression.close();
    }

    // Read lines from generated test output into a vector of strings
    std::vector<std::string> generatedLines;
    generatedLines.reserve(1535);
    {
      FileUtils::AutoGzIfstream fin_generated;
      fin_generated.openOrExit(params.outFileRoot + ".1.1.FastSMC.ibd.gz");
      for (std::string f_line; getline(fin_generated, f_line);) {
        generatedLines.emplace_back(f_line);
      }
      fin_generated.close();
    }

    REQUIRE(regressionLines.size() == generatedLines.size());

    for (auto lineNum = 0ul; lineNum < regressionLines.size(); ++lineNum) {
      REQUIRE(regressionLines.at(lineNum) == generatedLines.at(lineNum));
    }
  }
}
