#include <iostream>
#include <string>

#include "Data.hpp"
#include "DecodingParams.hpp"
#include "DecodingQuantities.hpp"
#include "FastSMC.hpp"
#include "HMM.hpp"
#include "Timer.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  srand(1234);

  // Parse input arguments
  DecodingParams params;
  if (!params.processCommandLineArgsFastSMC(argc, argv)) {
    cerr << "Error processing command line; exiting." << endl;
    exit(1);
  }

  DecodingQuantities decodingQuantities(params.decodingQuantFile.c_str());
  const auto sequenceLength = Data::countHapLines(params.inFileRoot);

  Data data(params.inFileRoot, sequenceLength, decodingQuantities.CSFSSamples, params.foldData, params.usingCSFS,
            params.jobInd, params.jobs);

  HMM hmm(data, decodingQuantities, params, !params.noBatches);
  hmm.decodeAll(params.jobs, params.jobInd);

  ASMC::RunFastSMC(params, data, hmm);
  return 0;
}
