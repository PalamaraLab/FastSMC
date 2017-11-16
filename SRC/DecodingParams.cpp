#include <iostream>
#include <string>
#include <sys/types.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include "DecodingParams.hpp"

using namespace std;

bool DecodingParams::processCommandLineArgs(int argc, char *argv[]) {

  namespace po = boost::program_options;

  string decodingModeString;

  po::options_description options;
  options.add_options()
  ("hapsFileRoot", po::value<string>(&hapsFileRoot)->required(),
   "prefix of hap|haps|hap.gz|haps.gz and sample|samples file")
  ("decodingQuantFile", po::value<string>(&decodingQuantFile),
   "decoding quantities file")
  ("jobs", po::value<int>(&jobs)->default_value(0),
   "number of jobs being done in parallel")
  ("jobInd", po::value<int>(&jobInd)->default_value(0),
   "job index (1..jobs)")
  ("outFileRoot", po::value<string>(&outFileRoot),
   "output file for sum of posterior distribution over pairs")
  ("mode", po::value<string>(&decodingModeString)->default_value("array"),
   "Decoding mode. Choose from {sequence, array}.")
  ("compress", po::bool_switch(&compress)->default_value(false),
   "Compress emission to binary (no CSFS)")
  // TODO: currently flipping major/minor when reading data. Instead, assume it's correctly coded to begin with
  ("useAncestral", po::bool_switch(&useAncestral)->default_value(false),
   "Assume ancestral alleles are coded as 1 in input (will assume 1 = minor otherwise)")
  // TODO: for debug. remove?
//  ("noBatches", po::bool_switch(&noBatches)->default_value(false),
//   "Do not decode in batches (for debugging, will remove)")
  ("skipCSFSdistance", po::value<float>(&skipCSFSdistance)->default_value(std::numeric_limits<float>::quiet_NaN()),
   "Genetic distance between two CSFS emissions")
  // TASKS
  ("posteriorSums", po::bool_switch(&doPosteriorSums)->default_value(false),
   "output file for sum of posterior distribution over pairs. O(sitesxstates) output")
  ("majorMinorPosteriorSums", po::bool_switch(&doMajorMinorPosteriorSums)->default_value(false),
    "output file for sum of posterior distribution over pairs, partitioned by major/minor alleles. O(sitesxstates) output")
  ("perPairMAP", po::bool_switch(&doPerPairMAP)->default_value(false),
   "output per-pair MAP at each site. O(sitesxsamples^2) output")
  ("perPairPosteriorMeans", po::value<string>(&expectedCoalTimesFile),
   "output per-pair posterior means at each site. Input: file containing expected coalescent times within each interval, one per line. O(sitesxsamples^2) output")
  ("withinOnly", po::bool_switch(&withinOnly)->default_value(false),
   "Only decode pairs within unphased individuals");

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
    cout << "skipCSFSdistance is " << skipCSFSdistance << endl;
    cout << "compress? " << compress << endl;
    cout << "useAncestral? " << useAncestral << endl;
    cout << "doPerPairPosteriorMean? " << doPerPairPosteriorMean << endl;
    cout << "doPerPairMAP? " << doPerPairMAP << endl;

    boost::algorithm::to_lower(decodingModeString);
    if (decodingModeString == string("sequence")) {
      decodingSequence = true;
      if (useAncestral) {
        decodingMode = DecodingMode::sequence;
        foldData = false;
      }
      else {
        decodingMode = DecodingMode::sequenceFolded;
        foldData = true;
      }
    } else if (decodingModeString == string("array")) {
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

    if (!doPosteriorSums && !doPerPairMAP && !doPerPairPosteriorMean && !doMajorMinorPosteriorSums) {
      cerr << "ERROR: At least one of --posteriorSums, --perPairMAP, --perPairPosteriorMean, --majorMinorPosteriorSums, must be specified"
           << endl;
      return false;
    }
  } catch (po::error &e) {
    cerr << "ERROR: " << e.what() << endl << endl;
    cerr << options << endl;
    return false;
  }

  return true;
}

