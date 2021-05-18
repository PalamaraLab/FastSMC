#include <iostream>
#include "Data.hpp"
#include "DecodingParams.hpp"
#include "DecodingQuantities.hpp"
#include "FileUtils.hpp"
#include "StringUtils.hpp"
#include "Timer.hpp"
#include "ASMC.hpp"
#include "HMM.hpp"

#include <fmt/format.h>


ASMC::ASMC::ASMC(DecodingParams params) : mParams{std::move(params)}, mData{mParams}, mHmm{mData, mParams}
{
}

ASMC::ASMC::ASMC(const std::string& inFileRoot, const std::string& decodingQuantFile, const std::string& outFileRoot)
    : mParams{inFileRoot,  decodingQuantFile,
              outFileRoot.empty() ? inFileRoot : outFileRoot, 1,
              1,           "array",
              false,       true,
              false,       false,
              0.f,         false,
              true,        false,
              "",          false,
              true,        true},
      mData{mParams}, mHmm{mData, mParams}
{
}

DecodingReturnValues ASMC::ASMC::decodeAllInJob() {
  std::cout << "Decoding job " << mParams.jobInd << " of " << mParams.jobs << "\n\n";

  std::cout << "Will decode " << mParams.decodingModeString << " data." << std::endl;
  std::cout << "Output will have prefix: " << mParams.outFileRoot << std::endl;

  if (mParams.compress) {
    std::cout << "Will use classic emission model (no CSFS)." << std::endl;
  } else {
    std::cout << "Minimum marker distance to use CSFS is set to " << mParams.skipCSFSdistance << "." << std::endl;
  }

  if (mParams.useAncestral) {
    std::cout << "Assuming ancestral alleles are correctly encoded." << std::endl;
  }

  if (mParams.doPosteriorSums) {
    std::cout << "Will output sum of posterior tables for all pairs." << std::endl;
  }

  if (mParams.doMajorMinorPosteriorSums) {
    std::cout << "Will output sum of posterior tables for all pairs, partitioned by major/minor alleles." << std::endl;
  }

  mHmm.decodeAll(mParams.jobs, mParams.jobInd);
  return mHmm.getDecodingReturnValues();
}

DecodingReturnValues ASMC::ASMC::decodePairs(const std::vector<uint>& individualsA,
                                             const std::vector<uint>& individualsB)
{
  mHmm.decodePairs(individualsA, individualsB);
  mHmm.finishDecoding();
  return mHmm.getDecodingReturnValues();
}
