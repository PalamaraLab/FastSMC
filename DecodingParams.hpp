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
    bool decodingSequence;
    bool foldData;
    bool usingCSFS;
    bool compress;
    bool useAncestral;
    float skipCSFSdistance;
    bool noBatches;

    // TASKS
    bool doPosteriorSums;
    bool doPerPairMAP; // output MAP for each pair
    bool doPerPairPosteriorMean; // output posterior mean for each pair
    string expectedCoalTimesFile; // expected coalescent times to output MAP/posterior mean times
    bool withinOnly; // only compute decoding within individuals
    //    float IBDthreshold;
    //    bool doMajorMinorPosteriorSums;
    //    bool doPerPairPosterior;
    //    bool doTrees;
    //    string indLabels;
    //    string secondHapsFileRoot;

    bool processCommandLineArgs(int argc, char *argv[]);

};

#endif
