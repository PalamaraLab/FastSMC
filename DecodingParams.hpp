#ifndef DECODINGPARAMS_HPP
#define DECODINGPARAMS_HPP

#include <string>

using namespace std;

enum class DecodingMode {
    sequenceFolded, arrayFolded, sequence, array
};

class DecodingParams {

public:

    string hapsFileRoot;
    string decodingQuantFile;
    string outFileRoot;
    int jobs;
    int jobInd;
    DecodingMode decodingMode;
    bool decodingSequence = false;
    bool foldData = false;
    bool usingCSFS = false;
    bool compress = false;
    bool useAncestral = false;
    float skipCSFSdistance;
    bool noBatches = false;

    // TASKS
    bool doPosteriorSums = false;
    bool doPerPairMAP = false; // output MAP for each pair
    bool doPerPairPosteriorMean = false; // output posterior mean for each pair
    string expectedCoalTimesFile; // expected coalescent times to output MAP/posterior mean times
    bool withinOnly = false; // only compute decoding within individuals
    //    float IBDthreshold;
    bool doMajorMinorPosteriorSums = false;
    //    bool doPerPairPosterior;
    //    bool doTrees;
    //    string indLabels;
    //    string secondHapsFileRoot;

    bool processCommandLineArgs(int argc, char *argv[]);

};

#endif
