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

#ifndef FASTSMC_DECODE_PAIRS_RETURN_STRUCT_HPP
#define FASTSMC_DECODE_PAIRS_RETURN_STRUCT_HPP

#include <Eigen/Core>
#include <fmt/format.h>

#include <vector>

struct DecodePairsReturnStruct {

private:
  std::size_t numWritten = 0ul;

  /// The pair indices: iInd, iHap, jInd, jHap
  Eigen::Array<int, Eigen::Dynamic, 4, Eigen::RowMajor> perPairIndices;

  /// The full set of posteriors: for each pair this is a (states * numSites) matrix
  Eigen::Array<Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Dynamic, 1> perPairPosteriors;

  /// The sum of all posteriors in perPairPosteriors: a (states * numSites) matrix
  Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> sumOfPosteriors;

  /// Posterior means: each row is an array of length numSites
  Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> perPairPosteriorMeans;

  Eigen::Array<float, 1, Eigen::Dynamic, Eigen::RowMajor> minPosteriorMeans;
  Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor> argminPosteriorMeans;

  Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> perPairMAPs;

  Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor> minMAPs;
  Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor> argminMAPs;

public:
  void resize(const std::vector<uint>& individualsA, const std::vector<uint>& individualsB, long int numSites,
              long int numStates)
  {
    Eigen::Index numPairsToDecode = 0ll;
    for (auto idx = 0ul; idx < individualsA.size(); ++idx) {
      if (individualsA[idx] == individualsB[idx]) {
        numPairsToDecode += 1ll;
      } else {
        numPairsToDecode += 4ll;
      }
    }

    perPairIndices.resize(numPairsToDecode, Eigen::NoChange);

    perPairPosteriors.resize(numPairsToDecode);
    for (Eigen::Index i = 0ll; i < perPairPosteriors.size(); ++i) {
      perPairPosteriors(i).resize(numStates, numSites);
    }

    sumOfPosteriors.resize(numStates, numSites);

    perPairPosteriorMeans.resize(numPairsToDecode, numSites);
    minPosteriorMeans.resize(numSites);
    argminPosteriorMeans.resize(numSites);

    perPairMAPs.resize(numPairsToDecode, numSites);
    minMAPs.resize(numSites);
    argminMAPs.resize(numSites);
  }

  void incrementNumWritten()
  {
    numWritten += 1;
  }

  void finaliseCalculations()
  {
    sumOfPosteriors.setZero();
    for (Eigen::Index pairIdx = 0ll; pairIdx < perPairPosteriors.size(); ++pairIdx) {
      sumOfPosteriors += perPairPosteriors(pairIdx);
    }

    for (Eigen::Index siteIdx = 0ll; siteIdx < perPairPosteriorMeans.cols(); ++siteIdx) {
      Eigen::Index argmin{};
      minPosteriorMeans(siteIdx) = perPairPosteriorMeans.col(siteIdx).minCoeff(&argmin);
      argminPosteriorMeans(siteIdx) = static_cast<int>(argmin);
    }

    for (Eigen::Index siteIdx = 0ll; siteIdx < perPairMAPs.cols(); ++siteIdx) {
      Eigen::Index argmin{};
      minMAPs(siteIdx) = perPairMAPs.col(siteIdx).minCoeff(&argmin);
      argminMAPs(siteIdx) = static_cast<int>(argmin);
    }
  }

  [[nodiscard]] const Eigen::Array<int, Eigen::Dynamic, 4, Eigen::RowMajor>& getPerPairIndices() const
  {
    return perPairIndices;
  }

  [[nodiscard]] const Eigen::Array<Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Dynamic, 1>& getPerPairPosteriors() const
  {
    return perPairPosteriors;
  }

  [[nodiscard]] const Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getSumOfPosteriors() const
  {
    return sumOfPosteriors;
  }

  [[nodiscard]] const Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getPerPairPosteriorMeans() const
  {
    return perPairPosteriorMeans;
  }

  [[nodiscard]] const Eigen::Array<float, 1, Eigen::Dynamic, Eigen::RowMajor>& getMinPosteriorMeans() const
  {
    return minPosteriorMeans;
  }

  [[nodiscard]] const Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor>& getArgminPosteriorMeans() const
  {
    return argminPosteriorMeans;
  }

  [[nodiscard]] const Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getPerPairMAPs() const
  {
    return perPairMAPs;
  }

  [[nodiscard]] const Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor>& getMinMAPs() const
  {
    return minMAPs;
  }

  [[nodiscard]] const Eigen::Array<int, 1, Eigen::Dynamic, Eigen::RowMajor>& getArgminMAPs() const
  {
    return argminMAPs;
  }

  [[nodiscard]] Eigen::Array<int, Eigen::Dynamic, 4, Eigen::RowMajor>& getModifiablePerPairIndices()
  {
    return perPairIndices;
  }

  [[nodiscard]] Eigen::Array<Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>, Eigen::Dynamic, 1>& getModifiablePerPairPosteriors()
  {
    return perPairPosteriors;
  }

  [[nodiscard]] Eigen::Array<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getModifiablePerPairPosteriorMeans()
  {
    return perPairPosteriorMeans;
  }

  [[nodiscard]] Eigen::Array<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>& getModifiablePerPairMAPs()
  {
    return perPairMAPs;
  }

  [[nodiscard]] std::size_t getNumWritten() const
  {
    return numWritten;
  }
};

#endif // FASTSMC_DECODE_PAIRS_RETURN_STRUCT_HPP
