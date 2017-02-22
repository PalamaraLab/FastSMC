#include <vector>
#include <string>
#include <cstdlib>
#include <iostream>

#include "Data.hpp"
#include "DecodingParams.hpp"
#include "DecodingQuantities.hpp"
#include "FileUtils.hpp"
// #include "Individual.hpp"
#include "StringUtils.hpp"
#include "Timer.hpp"
#include "HMM.cpp"

using namespace std;

int main(int argc, char *argv[]) {


    srand(1234);

    DecodingParams params;

    // parse input arguments
    if (!params.processCommandLineArgs(argc, argv)) {
        cerr << "Error processing command line; exiting." << endl;
        exit(1);
    }

    cout << "Output will have prefix: " << params.outFileRoot << endl;
    cout << "Output sum of posterior tables for all pairs? " << params.doPosteriorSums << endl;
    cout << "Use batches? " << !params.noBatches << endl;
    cout << "Only compute within individuals? " << params.withinOnly << endl;
    cout << "Output MAP at each site for all pairs? " << params.doPerPairMAP << endl;
    cout << "Output posterior mean at each site for all pairs? " << params.doPerPairPosteriorMean << endl;

    // used for benchmarking
    Timer timer;

    // read decoding quantities from file
    std::string str(params.decodingQuantFile.c_str());
    DecodingQuantities decodingQuantities(params.decodingQuantFile.c_str());
    printf("\n*** Read decoding quantities in %.3f seconds. ***\n\n", timer.update_time());
    cout << "CSFS samples: " << decodingQuantities.CSFSSamples << endl;

    int sequenceLength = Data::countHapLines(params.hapsFileRoot.c_str());
    Data data(params.hapsFileRoot.c_str(), sequenceLength, decodingQuantities.CSFSSamples, params.foldData, params.usingCSFS);
    printf("\n*** Read haps in %.3f seconds. ***\n\n", timer.update_time());

    // HMM hmm(data, decodingQuantities, DecodingParamsms.outFileRoot, params.doPerPairPosterior, params.doPosteriorSums, params.doMajorMinorPosteriorSums, params.doPerPairMAP, params.expectedCoalTimesFile, !params.doROHonly);
    HMM hmm(data, decodingQuantities, params, !params.noBatches);

    DecodingReturnValues decodingReturnValues = hmm.decodeAll(params.jobs, params.jobInd);
    vector < vector <float> > sumOverPairs = decodingReturnValues.sumOverPairs;

    // output sums over pairs (if requested)
    if (params.doPosteriorSums) {
        FileUtils::AutoGzOfstream fout; fout.openOrExit(params.outFileRoot + ".sumOverPairs.gz");
        for (int pos = 0; pos < data.sites; pos++) {
            for (uint k = 0; k < decodingQuantities.states; k++) {
                if (k) fout << "\t";
                fout << sumOverPairs[pos][k];
            }
            fout << endl;
        }
        fout.close();
    }
    if (params.doMajorMinorPosteriorSums) {
        vector < vector <float> > sumOverPairs00 = decodingReturnValues.sumOverPairs00;
        vector < vector <float> > sumOverPairs01 = decodingReturnValues.sumOverPairs01;
        vector < vector <float> > sumOverPairs11 = decodingReturnValues.sumOverPairs11;
        // Sum for 00
        FileUtils::AutoGzOfstream fout00;
        fout00.openOrExit(params.outFileRoot + ".00.sumOverPairs.gz");
        for (int pos = 0; pos < data.sites; pos++) {
            for (uint k = 0; k < decodingQuantities.states; k++) {
                if (k) fout00 << "\t";
                if (!data.siteWasFlippedDuringFolding[pos]) {
                    fout00 << sumOverPairs00[pos][k];
                }
                else {
                    fout00 << sumOverPairs11[pos][k];
                }
            }
            fout00 << endl;
        }
        fout00.close();
        // Sum for 01
        FileUtils::AutoGzOfstream fout01;
        fout01.openOrExit(params.outFileRoot + ".01.sumOverPairs.gz");
        for (int pos = 0; pos < data.sites; pos++) {
            for (uint k = 0; k < decodingQuantities.states; k++) {
                if (k) fout01 << "\t";
                fout01 << sumOverPairs01[pos][k];
            }
            fout01 << endl;
        }
        fout01.close();
        // Sum for 11
        FileUtils::AutoGzOfstream fout11;
        fout11.openOrExit(params.outFileRoot + ".11.sumOverPairs.gz");
        for (int pos = 0; pos < data.sites; pos++) {
            for (uint k = 0; k < decodingQuantities.states; k++) {
                if (k) fout11 << "\t";
                if (!data.siteWasFlippedDuringFolding[pos]) {
                    fout11 << sumOverPairs11[pos][k];
                }
                else {
                    fout11 << sumOverPairs00[pos][k];
                }
            }
            fout11 << endl;
        }
        fout11.close();
    }

    // cout << endl << "Upper-left 20 x 20 corner of sumOverPairs:" << endl;
    // for (uint i = 0; i < 20; i++) {
    //     for (uint j = 0; j < 20; j++)
    //         cout << sumOverPairs[i][j] << "\t";
    //     cout << endl;
    // }
}

