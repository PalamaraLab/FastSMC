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


#include <exception>
#include <iostream>
#include <string>
//#include <sys/types.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include "DecodingParams.hpp"
#include "Types.hpp"

using namespace std;

DecodingParams::DecodingParams()
      : hapsFileRoot("")
      , decodingQuantFile("")
      , outFileRoot(hapsFileRoot)
      , jobs(1)
      , jobInd(1)
      , decodingModeString("array")
      , decodingSequence(false)
      , usingCSFS(true)
      , compress(false)
      , useAncestral(false)
      , skipCSFSdistance(0.f)
      , noBatches(false)
      , doPosteriorSums(false)
      , doPerPairPosteriorMean(false)
      , expectedCoalTimesFile("")
      , withinOnly(false)
      , doMajorMinorPosteriorSums(false)
      , doPerPairMAP(false)
  {
  }

DecodingParams::DecodingParams(string _hapsFileRoot,
        string _decodingQuantFile,
        string _outFileRoot,
        int _jobs,
        int _jobInd,
        string _decodingModeString,
        bool _decodingSequence,
        bool _usingCSFS,
        bool _compress,
        bool _useAncestral,
        float _skipCSFSdistance,
        bool _noBatches,
        bool _doPosteriorSums,
        bool _doPerPairPosteriorMean,
        string _expectedCoalTimesFile,
        bool _withinOnly,
        bool _doMajorMinorPosteriorSums
        )
      : hapsFileRoot(_hapsFileRoot)
      , decodingQuantFile(_decodingQuantFile)
      , outFileRoot(_outFileRoot)
      , jobs(_jobs)
      , jobInd(_jobInd)
      , decodingModeString(_decodingModeString)
      , decodingSequence(_decodingSequence)
      , usingCSFS(_usingCSFS)
      , compress(_compress)
      , useAncestral(_useAncestral)
      , skipCSFSdistance(_skipCSFSdistance)
      , noBatches(_noBatches)
      , doPosteriorSums(_doPosteriorSums)
      , doPerPairPosteriorMean(_doPerPairPosteriorMean)
      , expectedCoalTimesFile(_expectedCoalTimesFile)
      , withinOnly(_withinOnly)
      , doMajorMinorPosteriorSums(_doMajorMinorPosteriorSums)
      , doPerPairMAP(false)
  {
     if(!processOptions()) throw std::exception();
  }


bool DecodingParams::processCommandLineArgs(int argc, char *argv[]) {

  namespace po = boost::program_options;

  po::options_description options;
  options.add_options()
  ("hapsFileRoot", po::value<string>(&hapsFileRoot)->required(),
   "Prefix of hap|haps|hap.gz|haps.gz and sample|samples file")
  ("decodingQuantFile", po::value<string>(&decodingQuantFile),
   "Decoding quantities file")
  ("jobs", po::value<int>(&jobs)->default_value(0),
   "Number of jobs being done in parallel")
  ("jobInd", po::value<int>(&jobInd)->default_value(0),
   "Job index (1..jobs)")
  ("outFileRoot", po::value<string>(&outFileRoot),
   "Output file for sum of posterior distribution over pairs (default: --hapsFileRoot argument)")
  ("mode", po::value<string>(&decodingModeString)->default_value("array"),
   "Decoding mode. Choose from {sequence, array}.")
  ("compress", po::bool_switch(&compress)->default_value(false),
   "Compress emission to binary (no CSFS)")
  // note: currently flipping major/minor when reading data. Instead, assume it's correctly coded to begin with
  ("useAncestral", po::bool_switch(&useAncestral)->default_value(false),
   "Assume ancestral alleles are coded as 1 in input (will assume 1 = minor otherwise)")
  // for debugging and pedagogical reasons.
 // ("noBatches", po::bool_switch(&noBatches)->default_value(false),
 //  "Do not decode in batches (for debugging, will remove)")
  ("skipCSFSdistance", po::value<float>(&skipCSFSdistance)->default_value(0.),
   "Genetic distance between two CSFS emissions")
  // main tasks
  ("majorMinorPosteriorSums", po::bool_switch(&doMajorMinorPosteriorSums)->default_value(false),
   "Output file for sum of posterior distribution over pairs, partitioned by major/minor alleles. O(sitesxstates) output")
  ("posteriorSums", po::bool_switch(&doPosteriorSums)->default_value(false),
   "Output file for sum of posterior distribution over pairs. O(sitesxstates) output");
  // ("perPairMAP", po::bool_switch(&doPerPairMAP)->default_value(false),
  //  "output per-pair MAP at each site. O(sitesxsamples^2) output")
  // ("perPairPosteriorMeans", po::value<string>(&expectedCoalTimesFile),
  //  "output per-pair posterior means at each site. O(sitesxsamples^2) output")
  // ("withinOnly", po::bool_switch(&withinOnly)->default_value(false),
  //  "Only decode pairs within unphased individuals");

  po::options_description visible("Options");
  visible.add(options);

  po::options_description all("All options");
  all.add(options);
  all.add_options()
  ("bad-args", po::value< vector <string> >(), "bad args")
  ;

  po::positional_options_description positional_desc;
  positional_desc.add("bad-args", -1); // for error-checking command line

  po::variables_map vm;
  po::command_line_parser cmd_line(argc, argv);
  cmd_line.options(all);
  cmd_line.style(po::command_line_style::default_style ^ po::command_line_style::allow_guessing);
  cmd_line.positional(positional_desc);

  try {
    po::store(cmd_line.run(), vm);

    po::notify(vm); // throws an error if there are any problems

    if (vm.count("bad-args")) {
      cerr << "ERROR: Unknown options:";
      vector <string> bad_args = vm["bad-args"].as< vector <string> >();
      for (uint i = 0; i < bad_args.size(); i++) cerr << " " << bad_args[i];
      cerr << endl;
      return false;
    }
  } catch (po::error &e) {
    cerr << "ERROR: " << e.what() << endl << endl;
    cerr << options << endl;
    return false;
  }

  if(processOptions()) {
    if (!doPosteriorSums && !doPerPairMAP && !doPerPairPosteriorMean && !doMajorMinorPosteriorSums) {
         cerr << "ERROR: At least one of --posteriorSums, --majorMinorPosteriorSums, must be specified"
              << endl;
         return false;
    }
  }
  else
      return false;
  return true;
}

bool DecodingParams::processOptions() {

    if (compress) {
      if (useAncestral) {
        cerr << "--compress & --useAncestral cannot be used together. A compressed emission cannot use ancestral allele information." << endl;
        exit(1);
      }
      if (!isnan(skipCSFSdistance)) {
        cerr << "--compress & --skipCSFSdistance cannot be used together. --compress is a shorthand for --skipCSFSdistance Infinity." << endl;
        exit(1);
      }
      skipCSFSdistance = std::numeric_limits<float>::infinity();
    }
    else {
      if (isnan(skipCSFSdistance)) {
        // default: use CSFS at all sites
        skipCSFSdistance = 0.f;
      }
    }

    if (!expectedCoalTimesFile.empty()) {
      doPerPairPosteriorMean = true;
    }

    if (skipCSFSdistance != std::numeric_limits<float>::infinity()) {
      usingCSFS = true;
    }

    boost::algorithm::to_lower(decodingModeString);
    if(decodingModeString == string("sequence"))
        decodingModeOverall = DecodingModeOverall::sequence;
    else if(decodingModeString == string("array"))
        decodingModeOverall = DecodingModeOverall::array;
    else {
       cerr << "Decoding mode should be one of {sequence, array}";
       return false;
    }

    if (decodingModeOverall == DecodingModeOverall::sequence) {
      decodingSequence = true;
      if (useAncestral) {
        decodingMode = DecodingMode::sequence;
        foldData = false;
      }
      else {
        decodingMode = DecodingMode::sequenceFolded;
        foldData = true;
      }
    } else if (decodingModeOverall == DecodingModeOverall::array) {
      decodingSequence = false;
      if (useAncestral) {
        decodingMode = DecodingMode::array;
        foldData = false;
      }
      else {
        decodingMode = DecodingMode::arrayFolded;
        foldData = true;
      }
    } else {
      cerr << "ERROR. Unknown decoding mode: " << decodingModeString << endl;
      exit(1);
    }

    if (decodingQuantFile.empty()) {
      cout << "Setting --decodingQuantFile to --hapsFileRoot + .decodingQuantities.bin" << endl;
      decodingQuantFile = hapsFileRoot + ".decodingQuantities.bin";
    }

    if ((jobs == 0) != (jobInd == 0)) {
      cerr << "ERROR: --jobs and --jobInd must either both be set or both be unset" << endl;
      return false;
    }

    if (jobs == 0) {
      jobs = 1;
      jobInd = 1;
    }

    if (jobInd <= 0 || jobInd > jobs) {
      cerr << "ERROR: --jobInd must be between 1 and --jobs inclusive" << endl;
      return false;
    }
    if (outFileRoot.empty()) {
      outFileRoot = hapsFileRoot;
      if (jobs > 0) {
        outFileRoot += "." + std::to_string(jobInd) + "-" + std::to_string(jobs);
      }
    }
   return true;
}

