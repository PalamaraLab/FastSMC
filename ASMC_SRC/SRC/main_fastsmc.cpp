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

vector<Individuals> all_ind;

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

int main(int argc, char* argv[])
{
  srand(1234);

  DecodingParams params;

  double TIME_start = get_cpu_time();
  double TIME_prev = TIME_start;
  float PAR_MIN_MAF = 0;
  float PAR_skip = 0;

  string line, discard;
  bool opt_error = 0;
  int c;

  // parse input arguments
  if (!params.processCommandLineArgsFastSMC(argc, argv)) {
    cerr << "Error processing command line; exiting." << endl;
    exit(1);
  }

  PAR_MIN_MATCH = params.min_m;
  PAR_MIN_MAF = params.min_maf;
  PAR_GAP = params.gap;
  PAR_skip = params.skip;
  MAX_seeds = params.max_seeds;

  // used for benchmarking
  Timer timer;
  // read decoding quantities from file
  std::string str(params.decodingQuantFile.c_str());
  DecodingQuantities decodingQuantities(params.decodingQuantFile.c_str());
  printf("\n*** Read decoding quantities in %.3f seconds. ***\n\n", timer.update_time());

  cout << "---------------------------" << endl;
  cout << "      	 READING DATA     	" << endl;
  cout << "---------------------------" << endl;
  cout << "CSFS samples: " << decodingQuantities.CSFSSamples << endl;
  unsigned int sequenceLength = Data::countHapLines(params.inFileRoot.c_str());

  // ifstream file_haps(params.inFileRoot + ".hap");
  FileUtils::AutoGzIfstream file_haps;
  if (FileUtils::fileExists(params.inFileRoot + ".hap.gz")) {
    file_haps.openOrExit(params.inFileRoot + ".hap.gz");
  } else if (FileUtils::fileExists(params.inFileRoot + ".hap")) {
    file_haps.openOrExit(params.inFileRoot + ".hap");
  } else if (FileUtils::fileExists(params.inFileRoot + ".haps.gz")) {
    file_haps.openOrExit(params.inFileRoot + ".haps.gz");
  } else if (FileUtils::fileExists(params.inFileRoot + ".haps")) {
    file_haps.openOrExit(params.inFileRoot + ".haps");
  } else {
    cerr << "ERROR. Could not find hap file in " + params.inFileRoot + ".hap.gz, " + params.inFileRoot + ".hap, " +
                params.inFileRoot + ".haps.gz, or " + params.inFileRoot + ".haps"
         << endl;
    exit(1);
  }

  // ifstream file_samp(params.inFileRoot + ".samples");
  FileUtils::AutoGzIfstream file_samp;
  if (FileUtils::fileExists(params.inFileRoot + ".samples")) {
    file_samp.openOrExit(params.inFileRoot + ".samples");
  } else if (FileUtils::fileExists(params.inFileRoot + ".sample")) {
    file_samp.openOrExit(params.inFileRoot + ".sample");
  } else {
    cerr << "ERROR. Could not find sample file in " + params.inFileRoot + ".sample or " + params.inFileRoot + ".samples"
         << endl;
    exit(1);
  }

  ifstream file_genm(params.map);
  if (!file_genm) {
    cerr << params.map << " could not be opened" << endl;
    return -1;
  }

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
    if (cur_g > 0 && (genetic_map[cur_g].first < genetic_map[cur_g - 1].first ||
                      genetic_map[cur_g].second < genetic_map[cur_g - 1].second)) {
      cerr << "ERROR: genetic map not in sorted order at line\n" << line << endl;
      return -1;
    }
    cur_g++;
  }
  file_genm.close();
  if (genetic_map.size() < 2) {
    cerr << "ERROR: genetic map must have at least two valid entries" << endl;
    return -1;
  }

  // cerr << "*** runtime : " << get_cpu_time() - TIME_start << "\t";
  cerr << genetic_map.size() << " genetic map entries read" << endl;
  printf("\n*** Read genetic map in %.3f seconds. ***\n\n", timer.update_time());

  Data data(params.inFileRoot, sequenceLength, decodingQuantities.CSFSSamples, params.foldData, params.usingCSFS,
                   params.jobInd, params.jobs, genetic_map);
  printf("\n*** Read haps in %.3f seconds. ***\n\n", timer.update_time());

  // freeing memory
  genetic_map.clear();
  vector<pair<unsigned long int, double>>().swap(genetic_map);

  cout << "---------------------------" << endl;
  cout << "   INFERRING IBD SEGMENTS  " << endl;
  cout << "---------------------------" << endl;

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
        all_ind.emplace_back(WORD_SIZE, CONST_READ_AHEAD, 2 * linectr);
        all_ind.emplace_back(WORD_SIZE, CONST_READ_AHEAD,2 * linectr + 1);
      } else {
        all_ind.emplace_back(WORD_SIZE, CONST_READ_AHEAD,2 * linectr);
        all_ind.emplace_back(WORD_SIZE, CONST_READ_AHEAD,2 * linectr);
      }
    }
    linectr++;
  }
  file_samp.close();
  num_ind_tot = data.sampleSize * 2;
  num_ind = all_ind.size();

  cerr << num_ind / 2 << " sample identifiers read" << endl;

  Marker cur_marker;
  // track position through genetic map
  cur_g = 0;
  unsigned int snp_ctr;
  char al[2], inp;
  string cur_al;

  // Storage for seeds
  SeedHash seeds;

  hash_size word[2];

  // Storage for extensions
  ExtendHash extend(WORD_SIZE, PAR_HAPLOID ? num_ind : num_ind / 2, PAR_HAPLOID);

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
  cerr << "processed " << GLOBAL_CURRENT_WORD * WORD_SIZE << " / " << all_markers->size() << " SNPs" << endl;

  printf("\n*** Inference done in %.3f seconds. ***\n\n", timer.update_time());

  return 0;
}
