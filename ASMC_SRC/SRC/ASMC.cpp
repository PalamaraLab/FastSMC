#include <iostream>
#include "Data.hpp"
#include "DecodingParams.hpp"
#include "DecodingQuantities.hpp"
#include "FileUtils.hpp"
#include "StringUtils.hpp"
#include "Timer.hpp"
#include "ASMC.hpp"
#include "HMM.hpp"


using namespace std;

DecodingReturnValues run(string haps_file_root, string decoding_quant_file,
         string out_file_root = "", DecodingModeOverall mode = DecodingModeOverall::array,
         int jobs = 0, int job_index = 0,
         float skip_csfs_distance = 0,
         bool compress = false, bool use_ancestral = false,
         bool posterior_sums = false, bool major_minor_posterior_sums = false) {

    srand(1234);

    DecodingParams params;
    params.hapsFileRoot = haps_file_root;
    params.decodingQuantFile = decoding_quant_file;
    params.outFileRoot = out_file_root.empty() ? haps_file_root : out_file_root;
    params.decodingModeOverall = mode;
    params.compress = compress;
    params.useAncestral = use_ancestral;
    params.skipCSFSdistance = skip_csfs_distance;
    params.doPosteriorSums = posterior_sums;
    params.doMajorMinorPosteriorSums = major_minor_posterior_sums;
    if(!params.processOptions()) {
        cerr << "Error in options processing" << endl;
        exit(1);
    }
    params.decodingModeString = params.decodingModeOverall == DecodingModeOverall::array ? "array" : "sequence";

    cout << "Decoding batch " << params.jobInd << " of " << params.jobs << "\n\n";

    cout << "Will decode " << params.decodingModeString << " data." << endl;
    cout << "Output will have prefix: " << params.outFileRoot << endl;
    if (params.compress)
        cout << "Will use classic emission model (no CSFS)." << endl;
    else
        cout << "Minimum marker distance to use CSFS is set to " << params.skipCSFSdistance << "." << endl;
    if (params.useAncestral)
        cout << "Assuming ancestral alleles are correctly encoded." << endl;
    if (params.doPosteriorSums)
        cout << "Will output sum of posterior tables for all pairs." << endl;
    if (params.doMajorMinorPosteriorSums)
        cout << "Will output sum of posterior tables for all pairs, partitioned by major/minor alleles." << endl;

    // if (params.noBatches)
    //     cout << "Will not process samples in batches (slower)." << endl;
    // if (!params.withinOnly)
    //     cout << "Will only decode maternal vs. paternal haplotypes." << endl;
    // if (params.doPerPairMAP)
    //     cout << "Will output MAP for all haploid pairs (DANGER: huge files)." << endl;
    // if (params.doPerPairPosteriorMean)
    //     cout << "Will output posterior mean for all haploid pairs (DANGER: huge files)." << endl;

    // used for benchmarking
    Timer timer;

    // read decoding quantities from file
    std::string str(params.decodingQuantFile.c_str());
    DecodingQuantities decodingQuantities(params.decodingQuantFile.c_str());
    printf("Read precomputed decoding info in %.9f seconds.\n", timer.update_time());
    // cout << "CSFS samples: " << decodingQuantities.CSFSSamples << endl;

    cout << "Data will be loaded from " << params.hapsFileRoot << "*\n";
    int sequenceLength = Data::countHapLines(params.hapsFileRoot.c_str());
    Data data(params.hapsFileRoot.c_str(), sequenceLength, decodingQuantities.CSFSSamples, params.foldData, params.usingCSFS);
    printf("Read haps in %.9f seconds.\n", timer.update_time());

    HMM hmm(data, decodingQuantities, params, !params.noBatches);

    hmm.decodeAll(params.jobs, params.jobInd);
    return hmm.getDecodingReturnValues();
}


