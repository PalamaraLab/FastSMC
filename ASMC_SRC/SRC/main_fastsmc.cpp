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

  const int WORD_SIZE = 64;
  const int CONST_READ_AHEAD = 10;
  const bool PAR_HAPLOID = true;
  ASMC::FastSMC fastSMC(WORD_SIZE, CONST_READ_AHEAD, PAR_HAPLOID);
  fastSMC.run(params, data, hmm);

  return 0;
}
