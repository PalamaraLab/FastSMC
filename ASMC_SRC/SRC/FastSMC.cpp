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

#include "FastSMC.hpp"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "FileUtils.hpp"
#include "HASHING/ExtendHash.hpp"
#include "HASHING/Individuals.hpp"
#include "HASHING/SeedHash.hpp"
#include "Timer.hpp"

ASMC::FastSMC::FastSMC(const int hashingWordSize, const int constReadAhead, const bool haploid)
    : mHashingWordSize{hashingWordSize}, mConstReadAhead{constReadAhead}, mHaploid{haploid}
{
}

void ASMC::FastSMC::run(const DecodingParams& params, const Data& data, HMM& hmm) {

  Timer timer;

  if(!params.GERMLINE) {
    std::cerr << "This method expects params.GERMLINE to be true. Exiting." << std::endl;
    exit(1);
  }

  std::vector<Individuals> all_ind;

  const auto jobs = params.jobs;
  const auto jobID = params.jobInd;
  const auto windowSize = data.windowSize;
  const auto w_i = data.w_i;
  const auto w_j = data.w_j;
  const bool is_j_above_diag = data.is_j_above_diag;

  auto isSampleInJob = [windowSize, w_i, w_j, jobs, jobID](unsigned i) {
    return ((i >= (uint)((w_i - 1) * windowSize) / 2 && i < (uint)(w_i * windowSize) / 2) ||
            (i >= (uint)((w_j - 1) * windowSize) / 2 && i < (uint)(w_j * windowSize) / 2) ||
            (jobs == jobID && i >= (uint)((w_j - 1) * windowSize) / 2));
  };

  // *** read Individuals information
  {
    FileUtils::AutoGzIfstream file_samp;
    file_samp.openOrExit(params.inFileRoot + ".samples");

    std::string line;
    std::stringstream ss;
    std::string map_field[2];

    unsigned int linectr = 0;
    while (getline(file_samp, line)) {
      std::vector<std::string> splitStr;
      std::istringstream iss(line);
      std::string buf;
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
        if (mHaploid) {
          all_ind.emplace_back(mHashingWordSize, mConstReadAhead, 2 * linectr);
          all_ind.emplace_back(mHashingWordSize, mConstReadAhead, 2 * linectr + 1);
        } else {
          all_ind.emplace_back(mHashingWordSize, mConstReadAhead,2 * linectr);
          all_ind.emplace_back(mHashingWordSize, mConstReadAhead,2 * linectr);
        }
      }

      linectr++;
    }
    file_samp.close();

    std::cerr << all_ind.size() / 2 << " sample identifiers read" << std::endl;

  }

  const auto PAR_MIN_MATCH = params.min_m;
  const auto PAR_MIN_MAF = params.min_maf;
  const auto PAR_GAP = params.gap;
  const auto PAR_skip = params.skip;
  const auto MAX_seeds = params.max_seeds;

  {
    FileUtils::AutoGzIfstream file_haps;
    file_haps.openOrExit(params.inFileRoot + ".hap.gz");

    const auto num_ind_tot = data.sampleSize * 2;
    const auto num_ind = all_ind.size();

    // Storage for seeds & extensions
    SeedHash seeds;
    ExtendHash extend(mHashingWordSize, num_ind, mHaploid);

    const std::vector<float>& all_markers = data.geneticPositions;

    int GLOBAL_READ_WORDS = 0;
    int GLOBAL_CURRENT_WORD = 0;

    std::string line;
    std::string map_field;
    std::string marker_id;
    std::stringstream ss;
    unsigned long marker_pos;

    // track position through genetic map
    unsigned int snp_ctr;
    char al[2], inp;

    while (true) {
      snp_ctr = 0;
      while (getline(file_haps, line)) {
        // read the meta data
        ss.clear();
        ss.str(line);
        ss >> map_field >> marker_id >> marker_pos >> al[0] >> al[1];
        if (map_field.empty()) {
          continue;
        }

        // restrict on MAF
        if (PAR_MIN_MAF > 0) {
          int maf_ctr = 0;
          for (int i = 0; i < num_ind_tot; i++) {
            ss >> inp;
            if (inp == '1') {
              maf_ctr++;
            }
          }
          auto maf = static_cast<float>(maf_ctr / static_cast<double>(num_ind_tot));
          if (maf < PAR_MIN_MAF || maf > 1 - PAR_MIN_MAF) {
            continue;
          }

          // re-load the data
          ss.clear();
          ss.str(line);
          ss >> map_field >> map_field >> map_field >> al[0] >> al[1];
        }

        // read haplotype
        unsigned int hap_ctr = 0;
        for (int i = 0; i < num_ind_tot; i++) {
          ss >> inp;
          if (isSampleInJob(i / 2)) {
            if (inp == '1') {
              all_ind[hap_ctr].setMarker(GLOBAL_READ_WORDS, snp_ctr);
            }
            hap_ctr++;
          }
        }
        snp_ctr++;

        if (snp_ctr % mHashingWordSize == 0) {
          if (++GLOBAL_READ_WORDS >= mConstReadAhead) {
            break;
          } else {
            std::cerr << "*** loading word buffer " << GLOBAL_READ_WORDS << " / " << mConstReadAhead << std::endl;
          }
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

      int GLOBAL_SKIPPED_WORDS = 0;
      int cur_seeds = seeds.size();

      // skip low-complexity words
      if (static_cast<float>(cur_seeds) / static_cast<float>(num_ind) > PAR_skip) {
        seeds.extendAllPairs(&extend, GLOBAL_CURRENT_WORD, all_ind, MAX_seeds, jobID, jobs, w_i, w_j, windowSize,
                             GLOBAL_READ_WORDS, GLOBAL_SKIPPED_WORDS, GLOBAL_CURRENT_WORD, is_j_above_diag);
        extend.clearPairsPriorTo(GLOBAL_CURRENT_WORD - PAR_GAP, GLOBAL_CURRENT_WORD, PAR_MIN_MATCH, all_markers, hmm);
      } else {
        std::cerr << "low complexity word - " << cur_seeds << " - skipping" << std::endl;
        extend.extendAllPairsTo(GLOBAL_CURRENT_WORD);
      }

      seeds.clear();

      for (unsigned int i = 0; i < num_ind; i++) {
        all_ind[i].clear(GLOBAL_CURRENT_WORD);
      }
      GLOBAL_CURRENT_WORD++;
    }

    extend.clearAllPairs(PAR_MIN_MATCH, all_markers, hmm);
    file_haps.close();

    hmm.finishFromGERMLINE();
    std::cerr << "processed " << GLOBAL_CURRENT_WORD * mHashingWordSize << " / " << all_markers.size() << " SNPs"
              << std::endl;
  }

  printf("\n*** Inference done in %.3f seconds. ***\n\n", timer.update_time());
}
