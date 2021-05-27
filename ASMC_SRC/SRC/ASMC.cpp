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

#include "ASMC.hpp"

#include <fmt/format.h>

#include <exception>
#include <iostream>

ASMC::ASMC::ASMC(DecodingParams params) : mParams{std::move(params)}, mData{mParams}, mHmm{mData, mParams}
{
}

ASMC::ASMC::ASMC(const std::string& inFileRoot, const std::string& decodingQuantFile, const std::string& outFileRoot)
    : mParams{inFileRoot,
              decodingQuantFile,
              outFileRoot.empty() ? inFileRoot : outFileRoot,
              1,
              1,
              "array",
              false,
              true,
              false,
              false,
              0.f,
              false,
              true,
              false,
              "",
              false,
              true,
              true},
      mData{mParams}, mHmm{mData, mParams}
{
}

DecodingReturnValues ASMC::ASMC::decodeAllInJob()
{
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

DecodePairsReturnStruct ASMC::ASMC::decodePairs(const std::vector<uint>& individualsA,
                                                const std::vector<uint>& individualsB)
{
  if (individualsA.empty() || individualsA.size() != individualsB.size()) {
    throw std::runtime_error(fmt::format("Expected two vectors of indices the same length, but got {} and {}",
                                         individualsA.size(), individualsB.size()));
  }

  mHmm.getDecodePairsReturnStruct().resize(individualsA, individualsB, mData.sites,
                                           mHmm.getDecodingQuantities().states);
  mHmm.setStorePerPairPosteriorMean(true);
  mHmm.setStorePerPairMap(true);
  mHmm.setStorePerPairPosterior(true);
  mHmm.decodePairs(individualsA, individualsB);
  mHmm.finishDecoding();
  mHmm.getDecodePairsReturnStruct().finaliseCalculations();
  return mHmm.getDecodePairsReturnStruct();
}
